#include "KinCat_ProcessDictionary.hpp"
#include <type_traits>

namespace KinCat {

template <typename DT>
ProcessDictionary<DT>::ProcessDictionary(const value_type_2d_view<ordinal_type, DT> &variant_orderings, //Lists the indecies of full configuration sites in symmetric orderings.
                                         const value_type_2d_view<ordinal_type, DT> &symmetry_orderings, //Lists the indecies of set subconfiguration sites in symmetric orderings. N.B. Indecies do NOT correspond to full configuration. Also note that set orderings are concatenated without any non-subconfiguration sites.   
                                         const value_type_2d_view<site_type, DT> &configurations,
                                         const value_type_1d_view<ordinal_type, DT> &configuration_sets,
                                         const value_type_2d_view<ordinal_type, DT> &processints,
                                         const value_type_2d_view<ordinal_type, DT> &constraints,
                                         const std::vector<std::vector<std::vector<ordinal_type>>> &process_constraints,
                                         const value_type_2d_view<real_type, DT> &rates,
                                         const value_type_2d_view<ordinal_type, DT> &process_symmetries)
    : _configurations(configurations), _configuration_sets(configuration_sets), _processints(processints),
      _constraints(constraints), 
      _rates(rates), _process_symmetries(process_symmetries) {
  //Read in specific sites for each configuration set. 
  ordinal_type n_sets = (_configuration_sets.extent(0) - 1);
  _config_sites = value_type_2d_view<ordinal_type, DT>(do_not_init_tag("config_sites"), n_sets, _configurations.extent(1)+1);
  ordinal_type conf_len(0);

  const auto config_sites_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _config_sites);
  const auto configurations_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _configurations);
  const auto configuration_sets_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _configuration_sets);

  for (ordinal_type i = 0, iend = n_sets; i < iend; ++i) {
    conf_len = 0;
    for (ordinal_type j = 0, jend = configurations_host.extent(1); j < jend; ++j) {
      config_sites_host(i,j+1) = 0; 
      if (configurations_host(configuration_sets_host(i), j) == -1) { 
        continue;
      }
      else {
        config_sites_host(i,conf_len+1) = j;
        conf_len++;
      }
    }
    config_sites_host(i,0) = conf_len;
  }
  Kokkos::deep_copy(_config_sites,config_sites_host);

  //_config_sites stores the indecies of the variant ordering identifying sites that are pertinent to the subconfiguration
  //The first entries (_config_sites(i,0)) stores the number of sites in the subconfiguration. 
  //_symmetry_sets stores the cumulative total of sites in the configuration sets.
  _symmetry_sets = value_type_1d_view<ordinal_type, DT>(do_not_init_tag("symmetry orderings sets indices"), n_sets);
  const auto symmetry_sets_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _symmetry_sets);

  symmetry_sets_host(0) = 0;
  for (ordinal_type i = 1, iend = n_sets; i < iend; ++i) {
    symmetry_sets_host(i) = symmetry_sets_host(i-1) + config_sites_host(i-1,0);
  }
  Kokkos::deep_copy(_symmetry_sets,symmetry_sets_host);

  value_type_2d_view<ordinal_type, DT> uniconfig_variants = variant_orderings; //uniconfig_variants is just a copy of the full configuration variants
  value_type_2d_view<ordinal_type, DT> sets_variants; //sets variants holds the configuration site indecies of the full configuration for each set. 
  sets_variants = value_type_2d_view<ordinal_type, DT>(do_not_init_tag("sets_variants"), uniconfig_variants.extent(0), n_sets*uniconfig_variants.extent(1));
  const auto sets_variants_host = Kokkos::create_mirror_view_and_copy(host_device_type(), sets_variants);
  const auto uniconfig_variants_host = Kokkos::create_mirror_view_and_copy(host_device_type(), uniconfig_variants);
  ordinal_type vidx(0), sidx(0), pidx(0), didx(0);
  std::vector<bool> site_bools(configurations_host.extent(1), true);

  const auto symmetry_orderings_host = Kokkos::create_mirror_view_and_copy(host_device_type(), symmetry_orderings);

  for (ordinal_type seti = 0; seti < n_sets; ++seti) {
    for (ordinal_type p = 0; p < uniconfig_variants.extent(0); ++p) {
      for (ordinal_type i = 0, iend = configurations_host.extent(1); i < iend; ++i) {
        site_bools[i] = true;
      }
      for (ordinal_type j = 0, jend = config_sites_host(seti, 0); j < jend; ++j) { //goes through set subconfiguration specific sites to assign orderings
        didx = symmetry_sets_host(seti) + j; //index of subconfiguration site in symmetry_orderings view. 
        pidx = config_sites_host(seti, symmetry_orderings_host(0, didx) + 1); // The variant ordering reference index of the subconfiguration site.
        sidx = config_sites_host(seti, symmetry_orderings_host(p, didx) + 1); // The variant ordering index of the subconfiguration site for the current variant
        //Since symmetry orderings of subconfigurations will replace full configuration variants, only use single variant ordering as stored in uniconfig_variants. 
        sets_variants_host(p, (seti * configurations_host.extent(1)) + pidx) = uniconfig_variants_host(0, sidx); //CHECK ROBUSTNESS IN OTHER SYSTEMS!!!
        site_bools[sidx] = false;
      }
      for (ordinal_type i = 0, iend = configurations_host.extent(1); i < iend; ++i) {
        if (site_bools[i]) {
          vidx = (seti * configurations_host.extent(1)) + i;
          sets_variants_host(p, vidx) = uniconfig_variants_host(0, i);
        }
      }
    }
  }
  _variant_orderings = value_type_2d_view<ordinal_type, DT> (
    do_not_init_tag("variant orderings"), sets_variants_host.extent(0), sets_variants_host.extent(1));
  const auto variant_orderings_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _variant_orderings);
  for (ordinal_type i = 0, iend = sets_variants_host.extent(0); i < iend; ++i) {
    for (ordinal_type j = 0, jend = sets_variants_host.extent(1); j < jend; ++j) {
      variant_orderings_host(i, j) = sets_variants_host(i,j);
    }
  }
  Kokkos::deep_copy(_variant_orderings, variant_orderings_host);

  _inverse_variant_orderings = value_type_2d_view<ordinal_type, DT>(
      do_not_init_tag("inverse variant orderings"), sets_variants_host.extent(0), sets_variants_host.extent(1));
  auto inverse_variant_orderings_host =
      Kokkos::create_mirror_view_and_copy(host_device_type(), _inverse_variant_orderings);
  //const auto variant_orderings_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _variant_orderings);
  for (ordinal_type i = 0, iend = variant_orderings_host.extent(0); i < iend; ++i) {
    for (ordinal_type seti = 0; seti < n_sets; ++seti) {
      for (ordinal_type j = 0, jend = (variant_orderings_host.extent(1)/n_sets); j < jend; ++j) {
        inverse_variant_orderings_host(i, (seti * jend) + variant_orderings_host(i, (seti * jend)+j)) = j;
      }
    }
  }
  Kokkos::deep_copy(_inverse_variant_orderings, inverse_variant_orderings_host);

  //This section initializes the config_first_ints view 
  _config_first_ints = value_type_1d_view<ordinal_type, DT>(
      do_not_init_tag("config_first_ints"), configurations_host.extent(0)+1);
  auto config_first_ints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _config_first_ints);
  auto processints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _processints);
  config_first_ints_host(0) = 0;
  ordinal_type current_config = 0;
  for (ordinal_type i = 1, iend = processints_host.extent(0); i < iend; ++i) {
    if (processints_host(i,0) ==  current_config) {
      continue;
    } else {
      while (processints_host(i,0) != current_config) {
        config_first_ints_host(current_config+1) = i;
        current_config++;
      }
    }
  }
  if (current_config <= configurations_host.extent(0)) {
    for (ordinal_type i = current_config+1, iend = configurations_host.extent(0)+1; i < iend; ++i) {
      config_first_ints_host(i) = processints_host.extent(0);
    }
  }
  Kokkos::deep_copy(_config_first_ints, config_first_ints_host);  
  
  _inverse_symmetry_orderings = value_type_2d_view<ordinal_type, DT>(
      do_not_init_tag("inverse symmetry orderings"), symmetry_orderings.extent(0), symmetry_orderings.extent(1));
  const auto inverse_symmetry_orderings_host = 
      Kokkos::create_mirror_view_and_copy(host_device_type(), _inverse_symmetry_orderings);
  for (ordinal_type i = 0, iend = symmetry_orderings_host.extent(0); i < iend; ++i) {
    for (ordinal_type j = 0, jend = symmetry_orderings_host.extent(1); j < jend; ++j) {
      inverse_symmetry_orderings_host(i, symmetry_orderings_host(i,j)) = j;
    }
  }
  Kokkos::deep_copy(_inverse_symmetry_orderings, inverse_symmetry_orderings_host);

  
  auto constraints_host_view = Kokkos::create_mirror_view_and_copy(host_device_type(), constraints);
  {
    for (ordinal_type eid = 0; eid < processints_host.extent(0); eid++) {
      const ordinal_type conf_a_idx = processints_host(eid, 0);
      const ordinal_type conf_b_idx = processints_host(eid, 1);
      const ordinal_type process_idx = processints_host(eid, 2);
      const ordinal_type n_variants = variant_orderings_host.extent(0);
      const ordinal_type n_sites_conf = variant_orderings_host.extent(1)/n_sets;
      ordinal_type seti = 0;
      for (ordinal_type i = 0; i < n_sets; ++i) {
        if (conf_a_idx < configuration_sets_host(i+1)) {
          seti = i; // conf_a_idx and conf_b_idx have the same set by definition
          break;
        }
      }

      const auto &constraints_test = process_constraints[process_idx];


      ordinal_type len_key = config_sites_host(seti, 0);

      value_type_1d_view<ordinal_type,DT> _config_sites_reverse(do_not_init_tag("_config_sites_reverse"), n_sites_conf);
      const auto config_sites_reverse_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _config_sites_reverse);
      for (ordinal_type i = 0; i < n_sites_conf; i++) {
        config_sites_reverse_host(i) = -1;
      }
      for (ordinal_type i = 0; i < len_key; i++) {
        config_sites_reverse_host(config_sites_host(seti, i + 1)) = i; 
      }

      using key_type = site_type[128];
      key_type a, a_variant, b, b_variant;
      for (ordinal_type k = 0; k < n_sites_conf; ++k) {
        a[k] = configurations_host(conf_a_idx, k);
        b[k] = configurations_host(conf_b_idx, k);
      }
      for (ordinal_type k = 1 ; k < n_sites_conf; ++k) {
      } 
      for (ordinal_type k = 1 ; k < n_sites_conf; ++k) {
      }

      for (ordinal_type i = 0; i < n_variants; ++i) {
        constraints_host_view(eid, i) = -1; /// initialize with an invalid value
        for (ordinal_type k = 0; k < len_key; ++k) {
          a_variant[k] = a[variant_orderings_host(i, config_sites_host(seti, k+1)+(seti*(config_sites_host.extent(1)-1)))];
        }

        /// check if input variant meets the constraints
        bool input_conf_variant_valid(true);
        for (ordinal_type l = 0, lend = constraints_test.size(); l < lend; ++l) {
          const auto &constraint = constraints_test[l];
          if (a_variant[config_sites_reverse_host(constraint[0])] != constraint[1]) {
            input_conf_variant_valid &= false;
          }
        }
        if (input_conf_variant_valid) {
          /// apply constraints to a_variant; b_variant matches to a_variant with constraints
          for (ordinal_type l = 0, lend = constraints_test.size(); l < lend; ++l) {
            const auto &constraint = constraints_test[l];
            a_variant[config_sites_reverse_host(constraint[0])] = constraint[2];
          }

          ordinal_type n_variants_satisfying_constraints(0);
          for (ordinal_type j = 0; j < n_variants; ++j) {
            for (ordinal_type k = 0; k < len_key; ++k) {
              b_variant[k] = b[variant_orderings_host(j, config_sites_host(seti, k+1)+(seti*(config_sites_host.extent(1)-1)))];
              }

            /// check if the final conf variant meets the constraints
            bool constraints_satisfied(true);
            for (ordinal_type k = 0; k < len_key; ++k)
              constraints_satisfied &= (b_variant[k] == a_variant[k]);

            if (constraints_satisfied) {
              if (n_variants_satisfying_constraints == 0)
                constraints_host_view(eid, i) = j;
              ++n_variants_satisfying_constraints;
            }
          }
        }
      }
    }
  }
  if (n_sets > 1)  {
    std::cout << "\n";
    std::cout << "Reset constraints (due to sets construction): " << constraints_host_view.extent(0) << ", " << constraints_host_view.extent(1) << "\n";
    for (ordinal_type i = 0, iend = constraints_host_view.extent(0); i < iend; ++i) {
      ordinal_type sum = (constraints_host_view(i, 0) < 0);
      std::cout << "  instance: " << i << ", " << processints_host(i, 2) << ", ["
                << constraints_host_view(i, 0);
      for (ordinal_type j = 1, jend = constraints_host_view.extent(1); j < jend; ++j) {
        std::cout << ", " << constraints_host_view(i, j);
        sum += (constraints_host_view(i, j) < 0);
      }
      std::cout << "]\n";
      if (sum == ordinal_type(constraints_host_view.extent(1)))
        std::cout << " <------------- this event is null event\n";
    }
  }
  Kokkos::deep_copy(_constraints, constraints_host_view);

  //Initialize timers
  //_functiontimer = value_type_1d_view<std::chrono::duration<double>, DT>("function timers", 10);

  //Initialize counters
  ordinal_type n_process_types = 0;
  for (ordinal_type i = 1, iend = processints_host.extent(0); i < iend; ++i) {
    if (processints_host(i,2) >= n_process_types) {
      n_process_types = processints_host(i,2) + 1;
    }
  }
  std::cout << "n_process_types found : " << n_process_types << '\n';
  _found_process_counter = value_type_1d_view<ordinal_type, DT> ("found process counters", n_process_types);
  _found_instance_counter = value_type_1d_view<ordinal_type, DT> ("found instance counters", constraints_host_view.extent(0));
}

template <typename DT>
void ProcessDictionary<DT>::overrideProcessRates(const std::vector<std::string> &process_override_array,
                                                 const std::vector<real_type> &process_rates_override_array,
                                                 const std::vector<ordinal_type> &instance_override_array,
                                                 const std::vector<real_type> &instance_rates_override_array,
                                                 const std::map<std::string, ordinal_type> &process_index,
                                                 const std::vector<std::string> &processes) {
  const ordinal_type n_samples = _rates.extent(0);
  const ordinal_type n_processints = _rates.extent(1);
  const ordinal_type n_processes = process_index.size();
  const ordinal_type n_processes_override = process_override_array.size();
  const ordinal_type n_processints_override = instance_override_array.size();
  std::vector<ordinal_type> process_index_indecies;

  std::vector<int>::iterator pit;
  for (auto& it : process_index) {
    if (process_index_indecies.size() == 0) {
      process_index_indecies.push_back(it.second);
      continue;
    }
    for (ordinal_type i = 0, iend = process_index_indecies.size(); i < iend; ++i) {
      if (it.second < process_index_indecies[i]) {
        pit = process_index_indecies.begin() + i;
        process_index_indecies.insert(pit, it.second);
        break;
      }
    }
    if (it.second > process_index_indecies[process_index_indecies.size()-1]) {
      process_index_indecies.push_back(it.second);
    }
  }

  //This section sets all samples in rates_override view with process specific rates from override array.
  value_type_2d_view<real_type, DT> rates_override(do_not_init_tag("rates override"), n_samples, n_processes);
  {
    const auto rates_override_host = Kokkos::create_mirror_view(host_device_type(), rates_override);
    Kokkos::deep_copy(rates_override_host, -1);
    for (ordinal_type j = 0, jend = n_processes_override; j < jend; ++j) {
      const ordinal_type i_jj = process_index.at(process_override_array[j]);
      ordinal_type jj = -1;
      for (ordinal_type i_idx = 0, i_idxend = n_processes; i_idx < i_idxend ; i_idx++) {
        if (i_jj == process_index_indecies[i_idx]) {
          jj = i_idx;
          break;
        }
      }
      
      if (jj == -1) {
        std::cout << "ERROR: Couldn't find process while overriding process rates!!! \n";
      }
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
    //This section writes any instance specific rates to override
    for (ordinal_type ij = 0 ; ij < rates.span(); ij++) {
          const ordinal_type i = ij / n_processints, j = ij % n_processints;
          auto process_name = processes[processints(j, 2)];
          ordinal_type process_type = process_index.at(process_name);
          for (ordinal_type i_idx = 0, i_idxend = n_processes; i_idx < i_idxend ; i_idx++) {
            if (process_type == process_index_indecies[i_idx]) {
              process_type = i_idx;
              break;
            }
          }
          real_type value = rates_override(i, process_type);
          for (ordinal_type k = 0, kend = n_processints_override; k < kend; ++k) {
            if (instance_override_array[k] == j) {
              value = instance_rates_override_array[k + (i * n_processints_override)];
              break;
            }
          }
          /// modified values are positive, unmodified values are -1
          /// overwrites rates view with the specified process and instance rates
          rates(i, j) = value < 0 ? rates(i, j) : value;
        }
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
  const auto configuration_sets_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _configuration_sets);
  const auto processints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _processints);
  const auto rates_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _rates);
  const auto constraints_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _constraints);
  
  ordinal_type n_errors(0);
  auto show_conf = [&](const ordinal_type i, const ordinal_type p) {
    if (p < 0) {
      std::cout << "[ " << i << ", [ " << ordinal_type(configurations_host(i, 0));
      for (ordinal_type j = 1, jend = configurations_host.extent(1); j < jend; ++j)
        std::cout << ", " << ordinal_type(configurations_host(i, j));
      std::cout << " ] ]";
    } else {
      std::cout << "[ " << i << ", [ " << ordinal_type(configurations_host(i, variant_orderings_host(p, 0)));
      for (ordinal_type j = 1, jend = configurations_host.extent(1); j < jend; ++j)
        std::cout << ", " << ordinal_type(configurations_host(i, variant_orderings_host(p, j)));
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
  /// 0. configuration sorted in ascending order (in unique sets)
  /// 1. configurations are unique
  {
    for (ordinal_type seti = 0; seti < (configuration_sets_host.extent(0)-1); ++seti) {
      for (ordinal_type i = configuration_sets_host(seti)+1, iend = configuration_sets_host(seti+1); i < iend; ++i) {
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
          ordinal_type seti;
          for (ordinal_type i_check = 1, iend_check = _config_first_ints.size(); i_check < iend_check; ++i_check) {
            if (idx_conf_a <  configuration_sets_host(i_check)) {
              seti = i_check-1;
              break;
            }
          }
          for (ordinal_type k = 0, kend = configurations_host.extent(1); k < kend; ++k) {
            atemp = variant_orderings_host(j, (seti * kend) + k);
            atemp = configurations_host(idx_conf_a, atemp);
            btemp = variant_orderings_host(constraint(j), (seti * kend) + k);
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
                  std::cout << "Constraint does not match\n";
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
    os << "-- inverse variant ordering : \n";
    const auto inverse_variant_orderings_host = Kokkos::create_mirror_view_and_copy(host_device_type(), _inverse_variant_orderings);
    for (ordinal_type i = 0, iend = inverse_variant_orderings_host.extent(0); i < iend; ++i) {
      os << indent4 << "[" << ordinal_type(inverse_variant_orderings_host(i, 0));
      for (ordinal_type j = 1, jend = inverse_variant_orderings_host.extent(1); j < jend; ++j) {
        os << ", " << inverse_variant_orderings_host(i, j);
      }
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
      os << "], [" << constraints_host(i, 0);
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
