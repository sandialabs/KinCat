#ifndef __KINCAT_PROCESS_DICTIONARY_HPP__
#define __KINCAT_PROCESS_DICTIONARY_HPP__

#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DeviceType> struct ProcessDictionary {
public:
  /// for GPU, we need to change this key type
  using key_type = site_type[128];

  /// configuration, process instances (these must be sorted)
  value_type_2d_view<ordinal_type, DeviceType> _variant_orderings;
  value_type_2d_view<ordinal_type, DeviceType> _inverse_variant_orderings;
  value_type_2d_view<ordinal_type, DeviceType> _symmetry_orderings;
  value_type_2d_view<ordinal_type, DeviceType> _inverse_symmetry_orderings;
  value_type_1d_view<ordinal_type, DeviceType> _symmetry_sets;
  value_type_2d_view<site_type, DeviceType> _configurations;
  value_type_1d_view<ordinal_type, DeviceType> _configuration_sets;
  value_type_1d_view<ordinal_type, DeviceType> _config_first_ints;
  value_type_2d_view<ordinal_type, DeviceType> _config_sites;
  value_type_2d_view<ordinal_type, DeviceType> _processints;
  value_type_2d_view<ordinal_type, DeviceType> _constraints;
  value_type_2d_view<real_type, DeviceType> _rates;
  value_type_2d_view<ordinal_type, DeviceType> _process_symmetries;

  /// timers
  value_type_1d_view<std::chrono::duration<double>, DeviceType> _functiontimer;
  value_type_1d_view<ordinal_type, DeviceType> _found_process_counter;
  value_type_1d_view<ordinal_type, DeviceType> _found_instance_counter;


public:
  KOKKOS_DEFAULTED_FUNCTION ProcessDictionary() = default;
  KOKKOS_DEFAULTED_FUNCTION ProcessDictionary(const ProcessDictionary &b) = default;
  ProcessDictionary(const value_type_2d_view<ordinal_type, DeviceType> &variant_orderings,
                    const value_type_2d_view<ordinal_type, DeviceType> &symmetry_orderings,
                    const value_type_2d_view<site_type, DeviceType> &configurations,
                    const value_type_1d_view<ordinal_type, DeviceType> &configuration_sets,
                    const value_type_2d_view<ordinal_type, DeviceType> &processints,
                    const value_type_2d_view<ordinal_type, DeviceType> &constraints,
                    const std::vector<std::vector<std::vector<ordinal_type>>> &process_constraints,
                    const value_type_2d_view<real_type, DeviceType> &rates,
                    const value_type_2d_view<ordinal_type, DeviceType> &process_symmetries);
  KOKKOS_DEFAULTED_FUNCTION ~ProcessDictionary() = default;

  void overrideProcessRates(const std::vector<std::string> &process_override_array,
                            const std::vector<real_type> &process_rates_override_array,
                            const std::vector<ordinal_type> &instance_override_array,
                            const std::vector<real_type> &instance_rates_override_array,
                            const std::map<std::string, ordinal_type> &process_index,
                            const std::vector<std::string> &processes);

  ordinal_type getNumberOfConfigurations() const;
  ordinal_type getNumberOfEvents() const;

  KOKKOS_INLINE_FUNCTION
  ordinal_type searchConfiguration(const ordinal_type seti, const ordinal_type p, const key_type &input_key) const{
    //This version searches for a single variant. Intended for sets version. Not currently used...
    ordinal_type idx_configuration = -1;
    ordinal_type c_count ;
    bool found_var(false);
        
        ordinal_type idx_config = 0;
        ordinal_type it = 0, step = 0;
        c_count = _configuration_sets(seti + 1)-_configuration_sets(seti);
        ordinal_type iter = 0;

        ordinal_type len_key = _config_sites(seti, 0);
        key_type key;
        for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
          key[i] = input_key[_variant_orderings(p, _config_sites(seti, i+1)+(seti*(_config_sites.extent(1)-1)))];
        }

        if (c_count >= 6) {
          //Binary ordered search through configuration list
          while (c_count > 0) {
            it = idx_config;
            step = c_count / 2;
            it += step;
            bool flag = false;
            for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
              const site_type a = _configurations((it+_configuration_sets(seti)), _config_sites(seti, i+1)), b = key[i];
              if (a == b)
                continue;
              flag = a < b;
              break;
            }
            if (flag) {
              idx_config = ++it;
              c_count -= step + 1;
            } else {
              c_count = step;
            }
            ++iter;
          }
        } else { //Brute force search through configuration list when short
          for (ordinal_type c_ind = 0; c_ind < c_count; c_ind++) {
            bool match = true;
            for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
              const site_type a = _configurations((c_ind+_configuration_sets(seti)), _config_sites(seti, i+1)), b = key[i];
              if (a != b) {
                match = false;
                break;
              }
            }
            if (match) {
              idx_config = c_ind;
              break;
            }
          }
        }

        found_var = true;
        for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
          const site_type a = _configurations((idx_config+_configuration_sets(seti)), _config_sites(seti, i+1)), b = key[i];
          found_var = (found_var && ((a == b) ));
        }
        if (found_var) {
          idx_configuration = idx_config+_configuration_sets(seti); // only define on first variant of configuration
        }
    return idx_configuration;
  }

  KOKKOS_INLINE_FUNCTION
  bool searchConfiguration(const key_type &input_key, ordinal_type &idx_configuration, 
                           ordinal_type &idx_variant) const {
    //This version is for single configuration set (i.e. UNICONFIG or FULLSYM). Faster since no view manipulation.
    bool found(false);
    for (ordinal_type p = 0, pend = _variant_orderings.extent(0); p < pend; ++p) {
      idx_configuration = 0;
      ordinal_type it = 0, step = 0, count = _configurations.extent(0), len_key = _configurations.extent(1);
      ordinal_type iter = 0;

      key_type key;
      for (ordinal_type i = 0, iend = len_key; i < iend; ++i)
        key[i] = input_key[_variant_orderings(p, i)];

      while (count > 0) {
        it = idx_configuration;
        step = count / 2;
        it += step;
        bool flag = false;
        for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
          const site_type a = _configurations(it, i), b = key[i];
          if (a == b)
            continue;
          flag = a < b;
          break;
        }
        if (flag) {
          idx_configuration = ++it;
          count -= step + 1;
        } else {
          count = step;
        }
        ++iter;
      }

      found = true;
      for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
        const site_type a = _configurations(idx_configuration, i), b = key[i];
        found = (found && (a == b));
      }
      if (found) {
        idx_variant = p;
        break;
      }
    }
    return found;
  }

  KOKKOS_INLINE_FUNCTION
  bool searchConfiguration(const key_type &input_key, 
                           Kokkos::Array<ordinal_type, 50> &idx_variant_array, Kokkos::Array<ordinal_type,50> & idx_configuration_array, const ordinal_type n_configs) const {
    
    bool all_found(false);
    ordinal_type configs_found = 0;
    
    for (ordinal_type seti = 0; seti < n_configs; seti++) {
      ordinal_type c_count ;
      for (ordinal_type p = 0, pend = _variant_orderings.extent(0); p < pend; ++p) {
        bool found_var(false);
        idx_variant_array[seti] = -1;

        
        ordinal_type idx_config = 0;
        ordinal_type it = 0, step = 0;
        c_count = _configuration_sets(seti + 1)-_configuration_sets(seti);
        ordinal_type iter = 0;

        ordinal_type len_key = _config_sites(seti, 0);
        key_type key;
        for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
          key[i] = input_key[_variant_orderings(p, _config_sites(seti, i+1)+(seti*(_config_sites.extent(1)-1)))];
        }

        if (c_count >= 6) {
          //Binary ordered search through configuration list
          while (c_count > 0) {
            it = idx_config;
            step = c_count / 2;
            it += step;
            bool flag = false;
            for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
              const site_type a = _configurations((it+_configuration_sets(seti)), _config_sites(seti, i+1)), b = key[i];
              if (a == b)
                continue;
              flag = a < b;
              break;
            }
            if (flag) {
              idx_config = ++it;
              c_count -= step + 1;
            } else {
              c_count = step;
            }
            ++iter;
          }
        } else { //Brute force search through configuration list when short
          for (ordinal_type c_ind = 0; c_ind < c_count; c_ind++) {
            bool match = true;
            for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
              const site_type a = _configurations((c_ind+_configuration_sets(seti)), _config_sites(seti, i+1)), b = key[i];
              if (a != b) {
                match = false;
                break;
              }
            }
            if (match) {
              idx_config = c_ind;
              break;
            }
          }
        }

        found_var = true;
        for (ordinal_type i = 0, iend = len_key; i < iend; ++i) {
          const site_type a = _configurations((idx_config+_configuration_sets(seti)), _config_sites(seti, i+1)), b = key[i];
          found_var = (found_var && ((a == b) ));
        }
        if (found_var) {
          idx_variant_array[seti] = p;
          idx_configuration_array[seti] = idx_config+_configuration_sets(seti); // only define on first variant of configuration 
          configs_found++; //This setup means that if error, configs_found will be index of first incorrectly found seti config.
          break;// break variants for loop
        }
      }
    }
    if (configs_found == n_configs) {
      all_found = true;
    }
    return all_found;
  }
  
  KOKKOS_INLINE_FUNCTION
  void update_timer( ordinal_type fid, std::chrono::duration<double> t_incr) const {
    _functiontimer(fid) += t_incr;
  }

  KOKKOS_INLINE_FUNCTION
  void update_process_counter( ordinal_type pid) const {
    _found_process_counter(pid) += 1;
  }

  KOKKOS_INLINE_FUNCTION
  void update_instance_counter( ordinal_type pid) const {
    _found_instance_counter(pid) += 1;
  }

  ordinal_type validateData(const std::vector<std::string> &processes,
                            const std::vector<std::vector<std::vector<ordinal_type>>> &process_constraints) const;
  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;
  std::ostream &showMe(std::ostream &os, const std::string &label, const std::vector<std::string> &processes,
                       const ordinal_type verbose = 0) const;
};

} // namespace KinCat

#endif
