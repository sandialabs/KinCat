#include "KinCat_Lattice.hpp"
#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DT>
Lattice<DT>::Lattice(const ordinal_type n_samples, const ordinal_type n_cells_x, const ordinal_type n_cells_y,
                     const value_type_2d_view<real_type, DT> edge_vectors,
                     const value_type_2d_view<real_type, DT> basis_vectors, const ordinal_type n_species)
    : _n_samples(n_samples), _n_cells_x(n_cells_x), _n_cells_y(n_cells_y), _n_basis(basis_vectors.extent(0)),
      _n_cells_domain_x(n_cells_x), _n_cells_domain_y(n_cells_y), _edge_vectors(edge_vectors),
      _basis_vectors(basis_vectors), _n_species(n_species), _sites(), _pattern() {
  _sites = value_type_2d_view<site_type, DT>("sites", _n_samples, _n_cells_x * _n_cells_y * _n_basis);
}

template <typename DT>
void Lattice<DT>::setDomain(const ordinal_type n_cells_domain_x, const ordinal_type n_cells_domain_y) {
  _n_cells_domain_x = n_cells_domain_x;
  _n_cells_domain_y = n_cells_domain_y;
}

template <typename DT>
void Lattice<DT>::createNeighborPattern(const value_type_2d_view<real_type, DT> &symmetry_operations,
                                        const value_type_2d_view<real_type, DT> &site_coordinates,
                                        const ordinal_type search_distance, const ordinal_type verbose) {
  const FunctionScope func_scope("Lattice::createNeighborPattern", __FILE__, __LINE__, verbose);

  Lattice<host_device_type> lattice_host;
  this->createMirrorObject(lattice_host);

  auto symmetry_operations_host = Kokkos::create_mirror_view_and_copy(host_device_type(), symmetry_operations);
  auto site_coordinates_host = Kokkos::create_mirror_view_and_copy(host_device_type(), site_coordinates);

  /// key types
  using coordinate_key_type = std::array<real_type, 2>;
  using pattern_key_type = std::array<real_type, 256>; /// array of coordinates
  KINCAT_CHECK_ERROR(site_coordinates.extent(0) > 128,
                     "Error: key array size is smaller than the number of site coordinates");

  /// compare function for floating point value key
  auto compare_coordinate_key = [](const coordinate_key_type &a, const coordinate_key_type &b) {
    const real_type eps(1.0e-8);
    bool r_val(false);
    for (ordinal_type i = 0; i < 2; ++i) {
      if (std::abs(a[i] - b[i]) <= eps)
        continue; /// they are same
      return (a[i] < b[i]);
    }
    return r_val;
  };
  auto compare_pattern_key = [](const pattern_key_type &a, const pattern_key_type &b) {
    const real_type eps(1.0e-8);
    bool r_val(false);
    for (ordinal_type i = 0; i < 256; ++i) {
      if (std::abs(a[i] - b[i]) <= eps)
        continue; /// they are same
      return (a[i] < b[i]);
    }
    return r_val;
  };

  /// create a map (key: supercell configuration, value: lattice location)
  std::map<coordinate_key_type, std::array<ordinal_type, 3>,
           std::function<bool(const coordinate_key_type &, const coordinate_key_type &)>>
      neighbor_coordinates_map(compare_coordinate_key);
  const ordinal_type dist = std::abs(search_distance);
  for (ordinal_type k0 = -dist; k0 <= dist; ++k0)
    for (ordinal_type k1 = -dist; k1 <= dist; ++k1)
      for (ordinal_type ib = 0; ib < _n_basis; ++ib) {
        real_type coord_x, coord_y;
        lattice_host.getCoordinates(k0, k1, ib, coord_x, coord_y);

        const coordinate_key_type key = {coord_x, coord_y};
        const std::array<ordinal_type, 3> val = {k0, k1, ib};
        neighbor_coordinates_map[key] = val;
      }

  if (verbose > 2) {
    std::cout << "  -- Neighbor coordinates (coordinates -> lattice index)\n";
    for (const auto &c : neighbor_coordinates_map) {
      const auto key = c.first;
      const auto val = c.second;
      std::cout << "    ";
      std::cout << "[" << key[0] << ", " << key[1] << "] -> "
                << " [ " << val[0] << ", " << val[1] << ", " << val[2] << "]\n";
    }
  }

  /// test for site coordinates (upto 64 site locations in a pattern)
  std::map<pattern_key_type, std::vector<ordinal_type>,
           std::function<bool(const pattern_key_type &, const pattern_key_type &)>>
      pattern_map(compare_pattern_key);

  /// rotate site coordinates
  for (ordinal_type j = 0, jend = symmetry_operations_host.extent(0); j < jend; ++j) {
    /// wrap the 1d array with 2x2 matrix
    value_type_2d_view<real_type, host_device_type> R(&symmetry_operations_host(j, 0), 2, 2);

    pattern_key_type key;
    for (ordinal_type i = 0, iend = site_coordinates_host.extent(0); i < iend; ++i) {
      const auto v = Kokkos::subview(site_coordinates_host, i, Kokkos::ALL());
      key[i * 2 + 0] = R(0, 0) * v(0) + R(0, 1) * v(1);
      key[i * 2 + 1] = R(1, 0) * v(0) + R(1, 1) * v(1);
    }
    pattern_map[key].push_back(j);
  }

  if (verbose > 2) {
    std::cout << "  -- Pattern (supercell coordinates)\n";
    for (const auto &c : pattern_map) {
      const auto key = c.first;
      const auto val = c.second;
      std::cout << "    key: [ [" << key[0] << ", " << key[1] << "]";
      for (ordinal_type i = 1, iend = site_coordinates_host.extent(0); i < iend; ++i)
        std::cout << ", [" << key[i * 2] << ", " << key[i * 2 + 1] << "]";
      std::cout << "], val: [" << val[0];
      for (ordinal_type i = 1, iend = val.size(); i < iend; ++i)
        std::cout << ", " << val[i];
      std::cout << "]\n";
    }
  }

  /// symmetry2pattern
  _symmetry2pattern =
      value_type_1d_view<ordinal_type, DT>(do_not_init_tag("symmetry pattern"), symmetry_operations_host.extent(0));
  auto symmetry2pattern_host = Kokkos::create_mirror_view(host_device_type(), _symmetry2pattern);
  Kokkos::deep_copy(symmetry2pattern_host, -1);
  {
    ordinal_type ip(0);
    for (const auto &c : pattern_map) {
      const auto val = c.second;
      for (ordinal_type i = 0, iend = val.size(); i < iend; ++i)
        symmetry2pattern_host(val[i]) = ip;
      ++ip;
    }

    for (ordinal_type i = 0, iend = symmetry2pattern_host.extent(0); i < iend; ++i) {
      KINCAT_CHECK_ERROR(symmetry2pattern_host(i) < 0,
                         "Error: symmetry2pattern initialization is failed (input is not consistent)");
    }
  }
  Kokkos::deep_copy(_symmetry2pattern, symmetry2pattern_host);

  if (verbose > 2) {
    std::cout << "  -- Symmetry2Pattern\n";
    std::cout << "    [" << symmetry2pattern_host(0);
    for (ordinal_type i = 1, iend = symmetry2pattern_host.extent(0); i < iend; ++i)
      std::cout << ", " << symmetry2pattern_host(i);
    std::cout << "]\n";
  }

  /// loop over the pattern map
  _pattern = value_type_3d_view<ordinal_type, DT>(do_not_init_tag("pattern"), pattern_map.size(),
                                                  site_coordinates.extent(0), 3);
  auto pattern_host = Kokkos::create_mirror_view(host_device_type(), _pattern);
  {
    const auto ev = lattice_host._edge_vectors;
    ordinal_type ip(0);
    for (const auto &p : pattern_map) {
      const auto p_key = p.first;
      const auto p_value = p.second;
      for (ordinal_type i = 0, iend = site_coordinates.extent(0); i < iend; ++i) {
        const real_type cx(p_key[i * 2 + 0]), cy(p_key[i * 2 + 1]);
        const coordinate_key_type n_key = {cx * ev(0, 0) + cy * ev(1, 0), cx * ev(0, 1) + cy * ev(1, 1)};
        const std::array<ordinal_type, 3> n_value = neighbor_coordinates_map[n_key];

        for (ordinal_type k = 0; k < 3; ++k)
          pattern_host(ip, i, k) = n_value[k];
      }
      ++ip;
    }
    if (verbose > 2) {
      std::cout << "  -- Pattern (lattice index)\n";
      for (ordinal_type k0 = 0, k0end = pattern_host.extent(0); k0 < k0end; ++k0) {
        std::cout << "    [ [ " << pattern_host(k0, 0, 0) << ", " << pattern_host(k0, 0, 1) << ", "
                  << pattern_host(k0, 0, 2) << "]";
        for (ordinal_type k1 = 1, k1end = pattern_host.extent(1); k1 < k1end; ++k1)
          std::cout << ", [" << pattern_host(k0, k1, 0) << ", " << pattern_host(k0, k1, 1) << ", "
                    << pattern_host(k0, k1, 2) << "]";
        std::cout << " ]\n";
      }
    }
  }
  Kokkos::deep_copy(_pattern, pattern_host);
}

template <typename DT> void Lattice<DT>::randomizeSites(const real_type fill, const ordinal_type seed) {
  const auto sites = _sites;

  Kokkos::deep_copy(sites, site_type(0));
  Kokkos::Random_XorShift64_Pool<DT> random(seed);

  const ordinal_type n_samples = sites.extent(0);
  const ordinal_type n_sites = sites.extent(1);
  const ordinal_type n_sites_fill = n_sites * fill;

  value_type_1d_view<ordinal_type, DT> site_locations(do_not_init_tag("site locations"), n_sites);
  {
    const auto site_locations_host = Kokkos::create_mirror_view(host_device_type(), site_locations);
    std::iota(site_locations_host.data(), site_locations_host.data() + n_sites, 0);

    std::mt19937 random_number_generator(seed);

    std::shuffle(site_locations_host.data(), site_locations_host.data() + n_sites, random_number_generator);
    Kokkos::deep_copy(site_locations, site_locations_host);
  }

  value_type_1d_view<site_type, DT> site_species(do_not_init_tag("site species"), n_sites);
  {
    const site_type aa = 1, bb = _n_species;
    Kokkos::fill_random(site_species, random, aa, bb); // site_type(1), site_type(_n_species));
  }

  using range_policy_type = Kokkos::RangePolicy<typename DT::execution_space>;
  range_policy_type policy(0, n_samples * n_sites);
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ordinal_type ij) {
        const ordinal_type i = ij / n_sites, j = ij % n_sites;
        if (j < n_sites_fill)
          sites(i, site_locations(j)) = site_species(j);
      });
}

template <typename DT> void Lattice<DT>::randomizeSites(const std::vector<real_type> &fill, const ordinal_type seed) {
  // Instantiates configuration based on random distribution of species in specific fill ratios.
  // Does not consider site permissions, could lead to configurations not allowed and not in dictionary. 
  const auto sites = _sites;
  const ordinal_type n_samples = sites.extent(0);
  const ordinal_type n_sites = sites.extent(1);
  const ordinal_type n_species = _n_species - 1;
    if (n_samples == 1) {
    KINCAT_CHECK_ERROR(fill.size() != (n_species), "Error: fill vector size does not match to the number of species");
  } else {
    KINCAT_CHECK_ERROR(fill.size() != (n_species * n_samples),
                       "Error: fill vector size does not match to the number of species and samples");
  }
  Kokkos::deep_copy(sites, site_type(0));
  Kokkos::Random_XorShift64_Pool<DT> random(seed);

  value_type_2d_view<ordinal_type, DT> site_locations(do_not_init_tag("site locations"), n_samples, n_sites);
  {
    std::mt19937 random_number_generator(seed);

    const auto site_locations_host = Kokkos::create_mirror_view(host_device_type(), site_locations);
    for (ordinal_type i = 0, iend = n_samples; i < iend; ++i) {
      const auto pbeg = site_locations_host.data() + i * n_sites;
      const auto pend = pbeg + n_sites;

      std::iota(pbeg, pend, 0);
      std::shuffle(pbeg, pend, random_number_generator);
    }
    Kokkos::deep_copy(site_locations, site_locations_host);
  }

  value_type_2d_view<real_type, DT> fill_ratio(do_not_init_tag("fill ratio"), n_samples, n_species);
  {
    const auto fill_ratio_host = Kokkos::create_mirror_view(host_device_type(), fill_ratio);
    real_type running_fill;
    for (ordinal_type i = 0, iend = fill_ratio_host.extent(0); i < iend; ++i) {
      running_fill = 0.0;
      for (ordinal_type j = 0; j < n_species; ++j) {
        fill_ratio_host(i, j) = fill[i * (n_species) + j];
        running_fill += fill_ratio_host(i, j);
      }
      KINCAT_CHECK_ERROR((running_fill > 1.0), "Sample fill ratios sum to more than 1.")
    }
    Kokkos::deep_copy(fill_ratio, fill_ratio_host);
  }

  using range_policy_type = Kokkos::RangePolicy<typename DT::execution_space>;
  range_policy_type policy(0, n_samples * n_sites);
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ordinal_type ij) {
        const ordinal_type i = ij / n_sites, j = ij % n_sites;
        ordinal_type site_start = 0;
        ordinal_type site_end = 0;
        for (int k = 1; k <= n_species; k++) {
          site_end = site_start + fill_ratio(i, k - 1) * n_sites;
          if ((j >= site_start) && (j < site_end)) {
            sites(i, site_locations(i, j)) = k;
          }
          site_start = site_end;
        }
      });
}

template <typename DT> void Lattice<DT>::readinSites(const value_type_2d_view<site_type, DT> &in_sites, const ordinal_type verbose) {
  //Written so that if single sample is read in, all samples are initialized with the same lattice.
  const FunctionScope func_scope("Lattice::readinSites", __FILE__, __LINE__, verbose);

  const auto sites = _sites;
  const ordinal_type n_samples = sites.extent(0);
  const ordinal_type n_sites = sites.extent(1);
  //Currently don't check if the species read-in are in the model. 
  const ordinal_type n_samples_in = in_sites.extent(1)/n_sites;

  KINCAT_CHECK_ERROR(((n_samples!=n_samples_in) && (n_samples_in!=1)), "Error: Number of read-in input samples does not match current ensemble.")

  Kokkos::deep_copy(sites, site_type(-1));
  const auto sites_host = Kokkos::create_mirror_view_and_copy(host_device_type(), sites);
  const auto in_sites_host = Kokkos::create_mirror_view_and_copy(host_device_type(), in_sites);

  if (n_samples_in == 1) {
    for (ordinal_type i = 0, iend = n_samples; i < iend; i++) {
      for (ordinal_type j = 0, jend = n_sites; j < jend; j++) {
        sites_host(i, j) = in_sites_host(0,j);
      } 
    }
  } else {
    for (ordinal_type i = 0, iend = n_samples; i < iend; i++) {
      for (ordinal_type j = 0, jend = n_sites; j < jend; j++) {
        sites_host(i, j) = in_sites_host(0,(i*n_sites)+j);
      }
    }
  }
  Kokkos::deep_copy(_sites, sites_host);
}

template <typename DT>
std::ostream &Lattice<DT>::showMe(std::ostream &os, const std::string &label, const ordinal_type verbose) const {
  Lattice<host_device_type> lattice_host;
  this->createMirrorObject(lattice_host);

  std::string indent2("  "), indent4("  "), indet6("      ");

  os << "-- Lattice : " << label << "\n";
  os << indent2 << "-- number of species : " << _n_species << "\n";
  os << indent2 << "-- domain : [" << _n_cells_x << ", " << _n_cells_y << ", " << _n_basis << "]\n";
  const auto ev = lattice_host._edge_vectors;
  if (ev.span()) {
    os << indent2 << "-- edge vectors : \n";
    os << indent4 << "[" << ev(0, 0) << "," << ev(0, 1) << "]\n";
    os << indent4 << "[" << ev(1, 0) << "," << ev(1, 1) << "]\n";
  }

  const auto bv = lattice_host._basis_vectors;
  if (bv.span()) {
    os << indent2 << "-- basis : \n";
    for (ordinal_type i = 0, iend = bv.extent(0); i < iend; ++i)
      os << indent4 << "[" << bv(i, 0) << "," << bv(i, 1) << "]\n";
  }

  const auto pattern = lattice_host._pattern;
  if (verbose > 0) {
    if (pattern.span()) {
      os << indent2 << "-- pattern : \n";
      for (ordinal_type k0 = 0, k0end = pattern.extent(0); k0 < k0end; ++k0) {
        os << indent4 << "[ [ " << pattern(k0, 0, 0) << ", " << pattern(k0, 0, 1) << ", " << pattern(k0, 0, 2) << "]";
        for (ordinal_type k1 = 1, k1end = pattern.extent(1); k1 < k1end; ++k1)
          os << ", [" << pattern(k0, k1, 0) << ", " << pattern(k0, k1, 1) << ", " << pattern(k0, k1, 2) << "]";
        os << " ]\n";
      }
    }
  }

  if (verbose > 1) {
    const auto s_host = lattice_host._sites;
    if (s_host.span()) {
      for (ordinal_type sid = 0; sid < _n_samples; ++sid) {
        const auto s = Kokkos::subview(s_host, sid, Kokkos::ALL());
        os << indent2 << "-- site sample : " << sid << "\n";
        for (ordinal_type i = 0, iend = s.extent(0); i < iend; ++i) {
          ordinal_type k0, k1, k2;
          real_type coord_x, coord_y;
          lattice_host.getLatticeIndex(i, k0, k1, k2);
          lattice_host.getCoordinates(k0, k1, k2, coord_x, coord_y);

          os << indent4 << "id : " << i << ", species : " << ordinal_type(s(i)) << ", coords : [" << coord_x << ","
             << coord_y << "]\n";
        }
      }
    }
  }
  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct Lattice<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct Lattice<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct Lattice<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat