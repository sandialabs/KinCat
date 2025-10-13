#ifndef __KINCAT_LATTICE_HPP__
#define __KINCAT_LATTICE_HPP__

#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DeviceType> struct Lattice {
public:
  /// sample
  ordinal_type _n_samples;

  /// domain information
  ordinal_type _n_cells_x, _n_cells_y, _n_basis;
  ordinal_type _n_cells_domain_x, _n_cells_domain_y;
  value_type_2d_view<real_type, DeviceType> _edge_vectors;
  value_type_2d_view<real_type, DeviceType> _basis_vectors;

  /// sites
  ordinal_type _n_species;
  value_type_2d_view<site_type, DeviceType> _sites;

  /// pattern
  value_type_3d_view<ordinal_type, DeviceType> _pattern;
  value_type_1d_view<ordinal_type, DeviceType> _symmetry2pattern;

  KOKKOS_DEFAULTED_FUNCTION Lattice() = default;
  KOKKOS_DEFAULTED_FUNCTION Lattice(const Lattice &b) = default;
  Lattice(const ordinal_type n_samples, const ordinal_type n_cells_x, const ordinal_type n_cells_y,
          const value_type_2d_view<real_type, DeviceType> edge_vectors,
          const value_type_2d_view<real_type, DeviceType> basis_vectors, const ordinal_type n_species);
  KOKKOS_DEFAULTED_FUNCTION ~Lattice() = default;
  void setDomain(const ordinal_type n_cells_domain_x, const ordinal_type n_cells_domain_y);

  void createNeighborPattern(const value_type_2d_view<real_type, DeviceType> &symmetry_operations,
                             const value_type_2d_view<real_type, DeviceType> &site_coordinates,
                             ordinal_type search_distance, ordinal_type verbose);
  void randomizeSites(const real_type fill, const ordinal_type seed);
  void randomizeSites(const std::vector<real_type> &fill, const ordinal_type seed);

  void readinSites(const value_type_2d_view<site_type, DeviceType> &in_sites,const ordinal_type verbose);

  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;

  template <typename ArgDeviceType> void copySites(const value_type_2d_view<site_type, ArgDeviceType> &sites) {
    Kokkos::deep_copy(_sites, sites);
  }
  template <typename ArgDeviceType>
  void copySites(const ordinal_type sid, const value_type_1d_view<site_type, ArgDeviceType> &sites) {
    Kokkos::deep_copy(Kokkos::subview(_sites, sid, Kokkos::ALL()), sites);
  }
  template <typename ArgDeviceType> void createMirrorObject(Lattice<ArgDeviceType> &mirror) const {
    mirror._n_cells_x = _n_cells_x;
    mirror._n_cells_y = _n_cells_y;
    mirror._n_basis = _n_basis;
    mirror._n_cells_domain_x = _n_cells_domain_x;
    mirror._n_cells_domain_y = _n_cells_domain_y;
    mirror._n_species = _n_species;

    using mirror_memory_space = typename ArgDeviceType::memory_space;
    mirror._edge_vectors = Kokkos::create_mirror_view_and_copy(mirror_memory_space(), _edge_vectors);
    mirror._basis_vectors = Kokkos::create_mirror_view_and_copy(mirror_memory_space(), _basis_vectors);
    mirror._sites = Kokkos::create_mirror_view_and_copy(mirror_memory_space(), _sites);
    mirror._pattern = Kokkos::create_mirror_view_and_copy(mirror_memory_space(), _pattern);
    mirror._symmetry2pattern = Kokkos::create_mirror_view_and_copy(mirror_memory_space(), _symmetry2pattern);
  }

  ///
  /// method to be used inside of a kernel
  ///
  KOKKOS_INLINE_FUNCTION
  ordinal_type getNumberOfSamples() const { return _n_samples; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type getNumberOfCells() const { return _n_cells_x * _n_cells_y; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type getNumberOfDomainCells() const { return _n_cells_domain_x * _n_cells_domain_y; }

  KOKKOS_INLINE_FUNCTION
  ordinal_type getNumberOfDomains(ordinal_type &n_domains_x, ordinal_type &n_domains_y) const {
    n_domains_x = _n_cells_x / _n_cells_domain_x;
    n_domains_y = _n_cells_y / _n_cells_domain_y;
    return n_domains_x * n_domains_y;
  }

  KOKKOS_INLINE_FUNCTION
  void getCellIndex(const ordinal_type k0, const ordinal_type k1, ordinal_type &cid) const {
    KINCAT_DEBUG_CHECK_ERROR(k0 < 0 || k0 >= _n_cells_x, "Error: k0 is out of bound, getCellIndex");
    KINCAT_DEBUG_CHECK_ERROR(k1 < 0 || k1 >= _n_cells_y, "Error: k1 is out of bound, getCellIndex");

    cid = k1 + k0 * (_n_cells_y);
  }

  KOKKOS_INLINE_FUNCTION
  void getDomainCellIndex(const ordinal_type k0, const ordinal_type k1, ordinal_type &d0, ordinal_type &d1,
                          ordinal_type &lid) const {
    KINCAT_DEBUG_CHECK_ERROR(k0 < 0 || k0 >= _n_cells_x, "Error: k0 is out of bound, getDomainCellIndex");
    KINCAT_DEBUG_CHECK_ERROR(k1 < 0 || k1 >= _n_cells_y, "Error: k1 is out of bound, getDomainCellIndex");

    d0 = k0 / _n_cells_domain_x;
    d1 = k1 / _n_cells_domain_y;

    const ordinal_type l0 = k0 - d0 * _n_cells_domain_x;
    const ordinal_type l1 = k1 - d1 * _n_cells_domain_y;

    lid = l1 + l0 * (_n_cells_domain_y);

  }

  KOKKOS_INLINE_FUNCTION
  void getLatticeCellIndex(const ordinal_type cid, ordinal_type &k0, ordinal_type &k1) const {
    k0 = cid / _n_cells_y;
    k1 = cid % _n_cells_y;

    KINCAT_DEBUG_CHECK_ERROR(k0 < 0 || k0 >= _n_cells_x, "Error: k0 is out of bound, getLatticeCellIndex");
    KINCAT_DEBUG_CHECK_ERROR(k1 < 0 || k1 >= _n_cells_y, "Error: k1 is out of bound, getLatticeCellIndex");
  }

  KOKKOS_INLINE_FUNCTION
  void getDomainLatticeCellIndex(const ordinal_type lid, ordinal_type &l0, ordinal_type &l1) const {
    l0 = lid / _n_cells_domain_y;
    l1 = lid % _n_cells_domain_y;

    KINCAT_DEBUG_CHECK_ERROR(l0 < 0 || l0 >= _n_cells_domain_x, "Error: l0 is out of bound");
    KINCAT_DEBUG_CHECK_ERROR(l1 < 0 || l1 >= _n_cells_domain_y, "Error: l1 is out of bound");
  }

  KOKKOS_INLINE_FUNCTION
  void getLatticeCellIndex(const ordinal_type d0, const ordinal_type d1, const ordinal_type l0, const ordinal_type l1,
                           ordinal_type &k0, ordinal_type &k1) const {
    k0 = d0 * _n_cells_domain_x + l0;
    k1 = d1 * _n_cells_domain_y + l1;

    KINCAT_DEBUG_CHECK_ERROR(k0 < 0 || k0 >= _n_cells_x, "Error: k0 is out of bound, getLatticeCellIndex");
    KINCAT_DEBUG_CHECK_ERROR(k1 < 0 || k1 >= _n_cells_y, "Error: k1 is out of bound, getLatticeCellIndex");
  }

  KOKKOS_INLINE_FUNCTION
  size_type getNumberOfSites() const { return _n_cells_x * _n_cells_y * _n_basis; }

  KOKKOS_INLINE_FUNCTION
  void getSiteIndex(const ordinal_type k0, const ordinal_type k1, const ordinal_type ib, ordinal_type &sid) const {
    KINCAT_DEBUG_CHECK_ERROR(k0 < 0 || k0 >= _n_cells_x, "Error: k0 is out of bound, getSiteIndex");
    KINCAT_DEBUG_CHECK_ERROR(k1 < 0 || k1 >= _n_cells_y, "Error: k1 is out of bound, getSiteIndex");
    KINCAT_DEBUG_CHECK_ERROR(ib < 0 || ib >= _n_basis, "Error: ib is out of bound");

    sid = ib + k1 * _n_basis + k0 * (_n_cells_y * _n_basis);
  }

  KOKKOS_INLINE_FUNCTION
  void getLatticeIndex(const ordinal_type sid, ordinal_type &k0, ordinal_type &k1, ordinal_type &ib) const {
    ordinal_type val(sid);
    {
      const ordinal_type stride = _n_basis * _n_cells_y;
      k0 = val / stride;
      val -= k0 * stride;
    }
    {
      const ordinal_type stride = _n_basis;
      k1 = val / stride;
      val -= k1 * stride;
    }
    ib = val;

    KINCAT_DEBUG_CHECK_ERROR(k0 < 0 || k0 >= _n_cells_x, "Error: k0 is out of bound, getLatticeIndex");
    KINCAT_DEBUG_CHECK_ERROR(k1 < 0 || k1 >= _n_cells_y, "Error: k1 is out of bound, getLatticeIndex");
    KINCAT_DEBUG_CHECK_ERROR(ib < 0 || ib >= _n_basis, "Error: ib is out of bound");
  }

  KOKKOS_INLINE_FUNCTION
  void adjustPeriodicBoundary(ordinal_type &k0, ordinal_type &k1) const {
    k0 = k0 < 0 ? _n_cells_x + k0 : k0;
    k0 = k0 >= _n_cells_x ? k0 - _n_cells_x : k0;

    k1 = k1 < 0 ? _n_cells_y + k1 : k1;
    k1 = k1 >= _n_cells_y ? k1 - _n_cells_y : k1;

    KINCAT_DEBUG_CHECK_ERROR(k0 < 0 || k0 >= _n_cells_x, "Error: k0 is out of bound, adjustPeriodicBoundary");
    KINCAT_DEBUG_CHECK_ERROR(k1 < 0 || k1 >= _n_cells_y, "Error: k1 is out of bound, adjustPeriodicBoundary");
  }

  KOKKOS_INLINE_FUNCTION
  void getCoordinates(const ordinal_type k0, const ordinal_type k1, const ordinal_type ib, real_type &coord_x,
                      real_type &coord_y) const {
    const auto ev = _edge_vectors;
    const auto bv = _basis_vectors;

    coord_x = (k0 + bv(ib, 0)) * ev(0, 0) + (k1 + bv(ib, 1)) * ev(1, 0);
    coord_y = (k0 + bv(ib, 0)) * ev(0, 1) + (k1 + bv(ib, 1)) * ev(1, 1);
  }
};

} // namespace KinCat

#endif