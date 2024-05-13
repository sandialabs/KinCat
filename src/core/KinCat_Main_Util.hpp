#ifndef __KINCAT_MAIN_UTIL_HPP__
#define __KINCAT_MAIN_UTIL_HPP__

#include "KinCat_Util.hpp"

#include "boost/json.hpp"
#include <regex>

namespace KinCat {

template <typename T>
ordinal_type parseArrayObject(const boost::json::object &obj, std::vector<ordinal_type> &shape, std::vector<T> &data) {
  shape = boost::json::value_to<std::vector<ordinal_type>>(obj.find("shape")->value());
  data = boost::json::value_to<std::vector<T>>(obj.find("data")->value());

  ordinal_type span(1);
  for (ordinal_type i = 0, iend = shape.size(); i < iend; ++i)
    span *= shape[i];

  KINCAT_CHECK_ERROR(span != data.size(), "Error: span given from shape information does not match to data span");
  return 0;
}

template <typename Ta, typename Tb>
ordinal_type copyStdVectorToKokkosView(const std::vector<Ta> &a, value_type_1d_view<Tb, host_device_type> &b) {
  if (a.size() > 0) {
    {
      std::stringstream ss;
      ss << "Error: dimension of a(" << a.size() << ") does not match to the dimension of b(" << b.extent(0) << ")\n";
      KINCAT_CHECK_ERROR(a.size() != b.extent(0), ss.str().c_str());
    }
    using range_policy_type = Kokkos::RangePolicy<typename host_device_type::execution_space>;
    range_policy_type policy(0, b.extent(0));
    Kokkos::parallel_for(policy, [&](const ordinal_type i) { b(i) = a[i]; });
  }
  return 0;
}

template <typename Ta, typename Tb>
ordinal_type copyStdVectorToKokkosView(const std::vector<std::vector<Ta>> &a,
                                       value_type_2d_view<Tb, host_device_type> &b) {
  if (a.size() > 0 && a[0].size() > 0) {
    {
      std::stringstream ss;
      ss << "Error: dimension of a(" << a.size() << ", " << a[0].size() << ") does not match to the dimension of b("
         << b.extent(0) << ", " << b.extent(1) << ")\n";
      KINCAT_CHECK_ERROR(a.size() != b.extent(0) && a[0].size() != b.extent(1), ss.str().c_str());
    }
    using range_policy_type = Kokkos::MDRangePolicy<typename host_device_type::execution_space, Kokkos::Rank<2>,
                                                    Kokkos::IndexType<ordinal_type>>;
    range_policy_type policy({0, 0}, {b.extent(0), b.extent(1)});
    Kokkos::parallel_for(policy, [&](const ordinal_type i, const ordinal_type j) { b(i, j) = a[i][j]; });
  }
  return 0;
}

template <typename Ta, typename Tb>
ordinal_type copyStdVectorToKokkosView(const std::vector<ordinal_type> &shape, const std::vector<Ta> &a,
                                       value_type_2d_view<Tb, host_device_type> &b) {
  if (a.size() > 0) {
    {
      KINCAT_CHECK_ERROR(shape.size() != 2, "Error: shape is not rank 2");

      std::stringstream ss;
      ss << "Error: dimension of a(" << shape[0] << ", " << shape[1] << ") does not match to the dimension of b("
         << b.extent(0) << ", " << b.extent(1) << ")\n";
      KINCAT_CHECK_ERROR(shape[0] != b.extent(0) && shape[1] != b.extent(1), ss.str().c_str());
    }
    using range_policy_type = Kokkos::MDRangePolicy<typename host_device_type::execution_space, Kokkos::Rank<2>,
                                                    Kokkos::IndexType<ordinal_type>>;
    range_policy_type policy({0, 0}, {b.extent(0), b.extent(1)});
    Kokkos::parallel_for(policy, [&](const ordinal_type i, const ordinal_type j) { b(i, j) = a[i * shape[1] + j]; });
  }
  return 0;
}

template <typename Ta, typename Tb>
ordinal_type copyStdVectorToKokkosView(const std::vector<ordinal_type> &shape, const std::vector<Ta> &a,
                                       value_type_3d_view<Tb, host_device_type> &b) {
  if (a.size() > 0) {
    {
      KINCAT_CHECK_ERROR(shape.size() != 3, "Error: shape is not rank 3");

      std::stringstream ss;
      ss << "Error: dimension of a(" << shape[0] << ", " << shape[1] << ", " << shape[2]
         << ") does not match to the dimension of b(" << b.extent(0) << ", " << b.extent(1) << ", " << b.extent(2)
         << ")\n";
      KINCAT_CHECK_ERROR(shape[0] != b.extent(0) && shape[1] != b.extent(1) && shape[2] != b.extent(2),
                         ss.str().c_str());
    }
    using range_policy_type = Kokkos::MDRangePolicy<typename host_device_type::execution_space, Kokkos::Rank<3>,
                                                    Kokkos::IndexType<ordinal_type>>;
    range_policy_type policy({0, 0, 0}, {b.extent(0), b.extent(1), b.extent(2)});
    Kokkos::parallel_for(policy, [&](const ordinal_type i, const ordinal_type j, const ordinal_type k) {
      b(i, j, k) = a[i * shape[1] * shape[2] + j * shape[2] + k];
    });
  }
  return 0;
}

template <typename T>
ordinal_type permute(const std::vector<ordinal_type> &index, const value_type_1d_view<T, host_device_type> &a) {
  value_type_1d_view<T, host_device_type> tmp(do_not_init_tag("tmp"), a.extent(0));
  Kokkos::deep_copy(tmp, a);
  using range_policy_type = Kokkos::RangePolicy<typename host_device_type::execution_space>;
  range_policy_type policy(0, tmp.extent(0));
  Kokkos::parallel_for(policy, [&](const ordinal_type i) { a(i) = tmp(index[i]); });
  return 0;
}
template <typename T>
ordinal_type permute(const std::vector<ordinal_type> &index, const value_type_2d_view<T, host_device_type> &a) {
  value_type_2d_view<T, host_device_type> tmp(do_not_init_tag("tmp"), a.extent(0), a.extent(1));
  Kokkos::deep_copy(tmp, a);
  using range_policy_type = Kokkos::RangePolicy<typename host_device_type::execution_space>;
  range_policy_type policy(0, tmp.extent(0));
  Kokkos::parallel_for(policy, [&](const ordinal_type i) {
    for (ordinal_type j = 0, jend = tmp.extent(1); j < jend; ++j)
      a(i, j) = tmp(index[i], j);
  });
  return 0;
}

std::string replaceEnvironmentVariable(const std::string &in);

ordinal_type parseInputfile(const std::string &filename, boost::json::value &jin, const ordinal_type verbose = 0);
ordinal_type validateDictionaryInput(const boost::json::value &tree, boost::json::object &root,
                                     const ordinal_type verbose = 0);
ordinal_type parseDictionaryInput(
    const boost::json::object &root, value_type_2d_view<real_type, device_type> &edge_vectors_view,
    value_type_2d_view<real_type, device_type> &basis_vectors_view,
    value_type_2d_view<real_type, device_type> &symmetry_operations_view,
    value_type_2d_view<real_type, device_type> &site_coordinates_view, ordinal_type &n_species,
    value_type_2d_view<ordinal_type, device_type> &variant_orderings_view, 
    value_type_2d_view<ordinal_type, device_type> &symmetry_orderings_view,
    ordinal_type &n_cells_interaction_x, ordinal_type &n_cells_interaction_y, 
    value_type_2d_view<site_type, device_type> &configurations_view,
    value_type_1d_view<ordinal_type, device_type> &configuration_sets,
    value_type_2d_view<ordinal_type, device_type> &instances_view,
    value_type_2d_view<ordinal_type, device_type> &constraints_view,
    value_type_2d_view<ordinal_type, device_type> &process_symmetries_view, std::vector<std::string> &processes,
    std::map<std::string, ordinal_type> &process_index,
    std::vector<std::vector<std::vector<ordinal_type>>> &process_constraints, const ordinal_type verbose = 0);
ordinal_type validateRatesInput(const boost::json::value &tree, boost::json::object &root,
                                const ordinal_type verbose = 0);
ordinal_type parseRatesInput(const boost::json::object &root, const std::vector<std::string> &processes,
                             const value_type_2d_view<ordinal_type, device_type> &instances_view,
                             value_type_2d_view<real_type, device_type> &rates_view, const ordinal_type verbose = 0);
ordinal_type sortAndRemapDictionary(const value_type_2d_view<site_type, device_type> &configurations_view,
                                    const value_type_1d_view<ordinal_type, device_type> &configuration_sets_view,
                                    const value_type_2d_view<ordinal_type, device_type> &instances_view,
                                    const value_type_2d_view<ordinal_type, device_type> &constraints_view,
                                    const value_type_2d_view<real_type, device_type> &rates_view,
                                    const ordinal_type verbose = 0);

ordinal_type validateSitesInput(const boost::json::value &tree, boost::json::object &root,
                                const ordinal_type verbose = 0);
ordinal_type parseSitesInput(const boost::json::object &root, std::string &site_init_type, ordinal_type &n_cells_x,
                             ordinal_type &n_cells_y, ordinal_type &n_basis, real_type &random_fill_ratio,
                             value_type_2d_view<site_type, device_type> &sites_view, value_type_1d_view< real_type, device_type> &t_sites,
                             const ordinal_type verbose = 0);
ordinal_type validateDumpInput(const boost::json::value &tree, boost::json::object &root,
                               const ordinal_type verbose = 0);

ordinal_type validateStatsInput(const boost::json::value &tree, boost::json::object &root,
                                const ordinal_type verbose = 0);

ordinal_type validateSolverInput(const boost::json::value &tree, boost::json::object &root,
                                 const ordinal_type verbose = 0);
ordinal_type parseSolverInput(const boost::json::object &root, std::string &solver_type, ordinal_type &n_cells_domain_x,
                              ordinal_type &n_cells_domain_y, ordinal_type &random_seed, ordinal_type &random_pool_size,
                              ordinal_type &n_kmc_kernel_launches, ordinal_type &n_kmc_steps_per_kernel_launch,
                              real_type &t_begin, real_type &t_end, real_type &t_dt, const ordinal_type verbose = 0);
ordinal_type validateEnsembleInput(const boost::json::value &tree, boost::json::object &root,
                                   const ordinal_type verbose = 0);
ordinal_type parseEnsembleInput(const boost::json::object &root, ordinal_type &n_samples,
                                bool &is_solver_random_number_variation, std::vector<real_type> &site_random_fill_ratio,
                                ordinal_type &rate_random_seed, std::vector<std::string> &process_rates_to_be_perturbed,
                                std::vector<real_type> &process_rates_scale_to_be_perturbed,
                                const ordinal_type verbose = 0);

struct Input {
public:
  /// lattice information
  value_type_2d_view<real_type, device_type> _edge_vectors;
  value_type_2d_view<real_type, device_type> _basis_vectors;
  value_type_2d_view<real_type, device_type> _symmetry_operations;
  value_type_2d_view<real_type, device_type> _site_coordinates;
  /// this include empty site as species zero
  ordinal_type _n_species;
  value_type_2d_view<ordinal_type, device_type> _variant_orderings;
  value_type_2d_view<ordinal_type, device_type> _symmetry_orderings;
  value_type_2d_view<site_type, device_type> _configurations;
  value_type_1d_view<ordinal_type, device_type> _configuration_sets;
  value_type_2d_view<ordinal_type, device_type> _processints;
  value_type_2d_view<ordinal_type, device_type> _constraints;
  value_type_2d_view<ordinal_type, device_type> _process_symmetries;
  std::vector<std::string> _processes;
  std::map<std::string, ordinal_type> _process_index;
  std::vector<std::vector<std::vector<ordinal_type>>> _process_constraints;

  /// rates
  value_type_2d_view<real_type, device_type> _rates;

  /// problem lattice
  ordinal_type _lattice_neighbor_search_space;
  ordinal_type _n_cells_x, _n_cells_y, _n_basis;

  /// site type: random, file
  std::string _site_init_type;

  /// random sites
  ordinal_type _site_random_seed;
  std::vector<real_type> _site_random_fill_ratio;

  /// sites from a file
  value_type_1d_view<real_type, device_type> _t_sites;
  value_type_2d_view<site_type, device_type> _sites;

  // dump
  std::string _dump_filename;
  std::string _restart_save_filename;
  real_type _dump_interval;

  /// statistics
  std::string _stats_filename;
  std::vector<std::string> _stats_list;
  real_type _stats_interval;

  /// kmc workflow
  std::string _solver_type;
  ordinal_type _solver_random_seed, _solver_random_pool_size;
  ordinal_type _n_kmc_kernel_launches, _n_kmc_steps_per_kernel_launch;
  ordinal_type _n_cells_domain_x, _n_cells_domain_y;
  ordinal_type _n_cells_interaction_x, _n_cells_interaction_y;

  real_type _t_begin, _t_end, _t_dt;

  /// ensemble
  ordinal_type _n_samples;
  std::vector<real_type> _site_random_fill_ratio_array;
  bool _is_solver_random_number_variation;
  std::vector<std::string> _process_override_array;
  std::vector<real_type> _process_rates_override_array;
  std::vector<ordinal_type> _instance_override_array;
  std::vector<real_type> _instance_rates_override_array;

  ordinal_type parse(const std::string input_filename, const ordinal_type verbose = 0);
};

} // namespace KinCat

#endif
