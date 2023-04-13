#include "KinCat_ProcessDictionary.hpp"
#include <type_traits>

namespace KinCat {

template <typename DT>
ProcessDictionary<DT>::ProcessDictionary(const value_type_2d_view<ordinal_type, DT> &variant_orderings,
                                         const value_type_2d_view<site_type, DT> &configurations,
                                         const value_type_2d_view<ordinal_type, DT> &processints,
                                         const value_type_2d_view<ordinal_type, DT> &constraints,
                                         const value_type_2d_view<real_type, DT> &rates,
                                         const value_type_2d_view<ordinal_type, DT> &process_symmetries)
    : _variant_orderings(variant_orderings), _configurations(configurations), _processints(processints),
      _constraints(constraints), _rates(rates), _process_symmetries(process_symmetries) {
  _inverse_variant_orderings = value_type_2d_view<ordinal_type, DT>(
      do_not_init_tag("inverse variant orderings"), variant_orderings.extent(0), variant_orderings.extent(1));
  const auto inverse_variant_orderings_host =
      Kokkos::create_mirror_view(host_device_type(), _inverse_variant_orderings);
  const auto variant_orderings_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _variant_orderings);
  for (ordinal_type i = 0, iend = variant_orderings_host.extent(0); i < iend; ++i)
    for (ordinal_type j = 0, jend = variant_orderings_host.extent(1); j < jend; ++j)
      inverse_variant_orderings_host(i, variant_orderings_host(i, j)) = j;
  Kokkos::deep_copy(_inverse_variant_orderings, inverse_variant_orderings_host);
}

template <typename DT>
void ProcessDictionary<DT>::overrideProcessRates(const std::vector<std::string> &process_override_array,
                                                 const std::vector<real_type> &process_rates_override_array,
                                                 const std::vector<ordinal_type> &instance_override_array,
                                                 const std::vector<real_type> &instance_rates_override_array,
                                                 const std::map<std::string, ordinal_type> &process_index) {
  const ordinal_type n_samples = _rates.extent(0);
  const ordinal_type n_processints = _rates.extent(1);
  const ordinal_type n_processes = process_index.size();
  const ordinal_type n_processes_override = process_override_array.size();
  const ordinal_type n_processints_override = instance_override_array.size();

  /// process specific variation
  value_type_2d_view<real_type, DT> rates_override(do_not_init_tag("rates override"), n_samples, n_processes);
  {
    const auto rates_override_host = Kokkos::create_mirror_view(host_device_type(), rates_override);
    Kokkos::deep_copy(rates_override_host, -1);
    for (ordinal_type j = 0, jend = n_processes_override; j < jend; ++j) {
      const ordinal_type jj = process_index.at(process_override_array[j]);
      for (ordinal_type i = 0, iend = n_samples; i < iend; ++i) {
        const real_type process_rate_at_ij = process_rates_override_array[j + (i * n_processes_override)];
        rates_override_host(i, jj) = process_rate_at_ij;
      }
    }
    Kokkos::deep_copy(rates_override, rates_override_host);
  }

  /// remap rates
  {
    const auto rates = _rates;
    const auto processints = _processints;

    using range_policy_type = Kokkos::RangePolicy<typename DT::execution_space>;
    range_policy_type policy(0, rates.span());
    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(const ordinal_type ij) {
          const ordinal_type i = ij / n_processints, j = ij % n_processints;
          real_type value = rates_override(i, processints(j, 2));
          for (ordinal_type k = 0, kend = n_processes_override; k < kend; ++k) {
            if (instance_override_array[k] == j) {
              value = instance_rates_override_array[j + (i * n_processints_override)];
              break;
            }
          }
          /// modified values are positive
          rates(i, j) = value < 0 ? rates(i, j) : value;
        });
    Kokkos::fence();
  }
}

template <typename DT> ordinal_type ProcessDictionary<DT>::getNumberOfConfigurations() const {
  return _configurations.extent(0);
}

template <typename DT> ordinal_type ProcessDictionary<DT>::getNumberOfEvents() const { return _processints.extent(0); }

template <typename DT>
ordinal_type ProcessDictionary<DT>::validateData(
    const std::vector<std::string> &processes,
    const std::vector<std::vector<std::vector<ordinal_type>>> &process_constraints) const {
  const auto variant_orderings_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _variant_orderings);
  const auto configurations_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _configurations);
  const auto processints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _processints);
  const auto rates_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _rates);
  const auto constraints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _constraints);

  ordinal_type n_errors(0);
  auto show_conf = [&](const ordinal_type i, const ordinal_type p) {
    if (p < 0) {
      std::cout << "[ " << i << ", [ " << configurations_host(i, 0);
      for (ordinal_type j = 1, jend = configurations_host.extent(1); j < jend; ++j)
        std::cout << ", " << configurations_host(i, j);
      std::cout << " ] ]";
    } else {
      std::cout << "[ " << i << ", [ " << configurations_host(i, variant_orderings_host(p, 0));
      for (ordinal_type j = 1, jend = configurations_host.extent(1); j < jend; ++j)
        std::cout << ", " << configurations_host(i, variant_orderings_host(p, j));
      std::cout << " ] ]";
    }
  };
  auto show_event = [&](const ordinal_type i) {
    {
      std::cout << "[ " << i << ", [ " << processints_host(i, 0) << ", " << processints_host(i, 1) << ", "
                << processints_host(i, 2) << "], ";
      std::cout << "process: " << processes[processints_host(i, 2)] << ", rate: " << rates_host(0, i);
      std::cout << " ] ]";
    }
  };

  /// validate configuration
  /// 0. configuration sorted in ascending order
  /// 1. configurations are unique
  {
    for (ordinal_type i = 1, iend = configurations_host.extent(0); i < iend; ++i) {
      const auto conf_a = Kokkos::subview(configurations_host, i - 1, Kokkos::ALL());
      const auto conf_b = Kokkos::subview(configurations_host, i, Kokkos::ALL());
      bool sorted_at_i = false;
      for (ordinal_type j = 0, jend = configurations_host.extent(1); j < jend; ++j) {
        if (conf_a(j) < conf_b(j)) {
          sorted_at_i = true;
          break;
        }
      }
      if (!sorted_at_i) {
        ++n_errors;
        std::cout << "Error: configurations are not sorted\n";
        show_conf(i - 1, -1);
        std::cout << ", ";
        show_conf(i, -1);
        std::cout << "\n";
      }
    }
  }

  /// validate process instances
  /// 0. instances are sorted in ascending order
  /// 1. rates should be positive (reverse process should be explicit)
  /// 2. an instance cannot map to the same initial configuration (null event)
  /// 3. variant map should satisfy the process constraints
  {
    for (ordinal_type i = 1, iend = processints_host.extent(0); i < iend; ++i) {
      const auto event_a = Kokkos::subview(processints_host, i - 1, Kokkos::ALL());
      const auto event_b = Kokkos::subview(processints_host, i, Kokkos::ALL());
      bool sorted_at_i = false;
      for (ordinal_type j = 0, jend = processints_host.extent(1); j < jend; ++j) {
        if (event_a(j) < event_b(j)) {
          sorted_at_i = true;
          break;
        }
      }
      if (!sorted_at_i) {
        ++n_errors;
        std::cout << "Error: process instances are not sorted\n";
        show_event(i - 1);
        std::cout << ", ";
        show_event(i);
        std::cout << "\n";
      }
    }
  }
  {
    for (ordinal_type i = 0, iend = processints_host.extent(0); i < iend; ++i) {
      const auto event = Kokkos::subview(processints_host, i, Kokkos::ALL());
      const auto constraint = Kokkos::subview(constraints_host, i, Kokkos::ALL());
      const auto rate = Kokkos::subview(rates_host, 0, i);

      const ordinal_type idx_conf_a = event(0);
      const ordinal_type idx_conf_b = event(1);
      const ordinal_type idx_process = event(2);

      if (rate() < 0) {
        ++n_errors;
        std::cout << "Error: process instance has a negative rate\n";
        show_event(i);
        std::cout << "\n";
      }

      {
        bool constraint_map_not_valid(true);
        for (ordinal_type j = 0, jend = constraint.extent(0); j < jend; ++j) {
          constraint_map_not_valid &= (constraint(j) < 0);
        }
        if (constraint_map_not_valid) {
          ++n_errors;
          std::cout << "Error: no valid variant map is available\n";
          show_event(i);
          std::cout << "\n";
        }
      }

      if (idx_conf_a == idx_conf_b) {
        /// same configuration; constraint should be applied and map to different varient
        bool constraint_map_valid(true);
        for (ordinal_type j = 0, jend = constraint.extent(0); j < jend; ++j) {
          constraint_map_valid &= (j != constraint(j));
        }
        if (!constraint_map_valid) {
          ++n_errors;
          std::cout << "Error: constraint maps to the same configuration and same variant\n";
          show_event(i);
          std::cout << "\n";
        }
      }

      for (ordinal_type j = 0, jend = constraint.extent(0); j < jend; ++j) {

        if (constraint(j) == -1) { // event not allowed from this configuration variant

        } else {
          // Check if can recover event constraints
          std::vector<std::vector<ordinal_type>> found_constraints;
          ordinal_type atemp, btemp;
          for (ordinal_type k = 0, kend = configurations_host.extent(1); k < kend; ++k) {
            atemp = variant_orderings_host(j, k);
            atemp = configurations_host(idx_conf_a, atemp);

            btemp = variant_orderings_host(constraint(j), k);
            btemp = configurations_host(idx_conf_b, btemp);

            if (atemp != btemp) { // Something changes at this site
              std::vector<ordinal_type> temp_constraint;
              temp_constraint.push_back(k);
              temp_constraint.push_back(atemp);
              temp_constraint.push_back(btemp);
              found_constraints.push_back(temp_constraint);
            }
          }
          bool error_flag = false;
          for (ordinal_type n = 0, nend = process_constraints[idx_process].size(); n < nend; ++n) {
            for (ordinal_type k = 0, kend = found_constraints.size(); k < kend; ++k) {
              if (process_constraints[idx_process][n][0] == found_constraints[k][0]) {
                if ((process_constraints[idx_process][n][1] == found_constraints[k][1]) &&
                    (process_constraints[idx_process][n][2] == found_constraints[k][2])) {
                } else {
                  std::cout << "Constraint does not match";
                  error_flag = true;
                }
                break;
              }
            }
            // Not in found constraints, check if bystander
            if (error_flag == true &&
                (process_constraints[idx_process][n][1] != process_constraints[idx_process][n][2])) {
              std::cout << "Not a bystander\n";
              error_flag = true; // Not a bystander constraint
            }
          }
          if (error_flag) {
            ++n_errors;
            std::cout << "Error: process instance does not match input constraints ";
            show_event(i);
            std::cout << "\n";
          }
        }
      }
    }
  }
  return n_errors;
}

template <typename DT>
std::ostream &ProcessDictionary<DT>::showMe(std::ostream &os, const std::string &label,
                                            const ordinal_type verbose) const {
  const std::string indent2("  "), indent4("  "), indet6("      ");
  os << "-- ProcessDictionary : " << label << "\n";
  os << indent2 << "-- variant ordering : (" << _variant_orderings.extent(0) << ", " << _variant_orderings.extent(1)
     << ")\n";
  if (verbose > 0) {
    const auto variant_orderings_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _variant_orderings);
    for (ordinal_type i = 0, iend = variant_orderings_host.extent(0); i < iend; ++i) {
      os << indent4 << "[" << ordinal_type(variant_orderings_host(i, 0));
      for (ordinal_type j = 1, jend = variant_orderings_host.extent(1); j < jend; ++j)
        os << ", " << variant_orderings_host(i, j);
      os << "]\n";
    }
  }

  os << indent2 << "-- configuration list : (" << _configurations.extent(0) << ", " << _configurations.extent(1)
     << ")\n";
  if (verbose > 0) {
    const auto configurations_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _configurations);
    for (ordinal_type i = 0, iend = configurations_host.extent(0); i < iend; ++i) {
      os << indent4 << "[" << i << ", [" << ordinal_type(configurations_host(i, 0));
      for (ordinal_type j = 1, jend = configurations_host.extent(1); j < jend; ++j)
        os << ", " << ordinal_type(configurations_host(i, j));
      os << "] ]\n";
    }
  }

  os << indent2 << "-- process instance list : (" << _processints.extent(0) << ", " << _processints.extent(1)
     << "), constraints: (" << _constraints.extent(0) << ", " << _constraints.extent(1) << "), rates : ("
     << _rates.extent(1) << ")\n";
  if (verbose > 0) {
    const auto processints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _processints);
    const auto rates_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _rates);
    const auto constraints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _constraints);
    for (ordinal_type i = 0, iend = processints_host.extent(0); i < iend; ++i) {
      os << indent4 << "[" << i << ", [" << processints_host(i, 0);
      for (ordinal_type j = 1, jend = processints_host.extent(1); j < jend; ++j) {
        os << ", " << processints_host(i, j);
      }
      os << "], [" << _constraints(i, 0);
      for (ordinal_type j = 1, jend = constraints_host.extent(1); j < jend; ++j) {
        os << ", " << constraints_host(i, j);
      }
      os << "], " << rates_host(0, i) << "]\n";
    }
  }
  return os;
}

template <typename DT>
std::ostream &ProcessDictionary<DT>::showMe(std::ostream &os, const std::string &label,
                                            const std::vector<std::string> &processes,
                                            const ordinal_type verbose) const {
  const std::string indent2("  "), indent4("  "), indet6("      ");
  os << "-- ProcessDictionary : " << label << "\n";

  os << indent2 << "-- process instance list : (" << _processints.extent(0) << "), process list : (" << processes.size()
     << ")\n";
  if (verbose > 0) {
    const auto configurations_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _configurations);
    const auto processints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _processints);
    const auto rates_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _rates);

    for (ordinal_type i = 0, iend = processints_host.extent(0); i < iend; ++i) {
      const ordinal_type idx_conf_a = processints_host(i, 0);
      const ordinal_type idx_conf_b = processints_host(i, 1);
      const ordinal_type idx_process = processints_host(i, 2);

      os << indent4 << "[" << i << ", [ [" << idx_conf_a << ", [" << ordinal_type(configurations_host(idx_conf_a, 0));
      for (ordinal_type j = 1, jend = configurations_host.extent(1); j < jend; ++j)
        os << ", " << ordinal_type(configurations_host(idx_conf_a, j));
      os << "] ],  [";
      os << idx_conf_b << ", [" << ordinal_type(configurations_host(idx_conf_b, 0));
      for (ordinal_type j = 1, jend = configurations_host.extent(1); j < jend; ++j)
        os << ", " << ordinal_type(configurations_host(idx_conf_b, j));
      os << "] ] ], " << rates_host(0, i) << ", [" << idx_process << ", \"" << processes[idx_process] << "\"] ]\n";
    }
  }
  return os;
}

/// eti for individual device type
#if defined(KOKKOS_ENABLE_SERIAL)
template struct ProcessDictionary<typename UseThisDevice<Kokkos::Serial>::type>;
#endif
#if defined(KOKKOS_ENABLE_OPENMP)
template struct ProcessDictionary<typename UseThisDevice<Kokkos::OpenMP>::type>;
#endif
#if defined(KOKKOS_ENABLE_CUDA)
template struct ProcessDictionary<typename UseThisDevice<Kokkos::Cuda>::type>;
#endif

} // namespace KinCat
