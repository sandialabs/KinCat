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
  value_type_2d_view<site_type, DeviceType> _configurations;
  value_type_2d_view<ordinal_type, DeviceType> _processints;
  value_type_2d_view<ordinal_type, DeviceType> _constraints;
  value_type_2d_view<real_type, DeviceType> _rates;
  value_type_2d_view<ordinal_type, DeviceType> _process_symmetries;

public:
  KOKKOS_DEFAULTED_FUNCTION ProcessDictionary() = default;
  KOKKOS_DEFAULTED_FUNCTION ProcessDictionary(const ProcessDictionary &b) = default;
  ProcessDictionary(const value_type_2d_view<ordinal_type, DeviceType> &variant_orderings,
                    const value_type_2d_view<site_type, DeviceType> &configurations,
                    const value_type_2d_view<ordinal_type, DeviceType> &processints,
                    const value_type_2d_view<ordinal_type, DeviceType> &constraints,
                    const value_type_2d_view<real_type, DeviceType> &rates,
                    const value_type_2d_view<ordinal_type, DeviceType> &process_symmetries);
  KOKKOS_DEFAULTED_FUNCTION ~ProcessDictionary() = default;

  void overrideProcessRates(const std::vector<std::string> &process_override_array,
                            const std::vector<real_type> &process_rates_override_array,
                            const std::vector<ordinal_type> &instance_override_array,
                            const std::vector<real_type> &instance_rates_override_array,
                            const std::map<std::string, ordinal_type> &process_index);

  ordinal_type getNumberOfConfigurations() const;
  ordinal_type getNumberOfEvents() const;

  KOKKOS_INLINE_FUNCTION
  bool searchConfiguration(const key_type &input_key, ordinal_type &idx_configuration,
                           ordinal_type &idx_variant) const {
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
  bool searchFirstEvent(const ordinal_type idx_configuration, ordinal_type &idx_first_event) const {
    ordinal_type it = 0, step = 0, count = _processints.extent(0);
    while (count > 0) {
      it = idx_first_event;
      step = count / 2;
      it += step;
      if (_processints(it, 0) < idx_configuration) {
        idx_first_event = ++it;
        count -= step + 1;
      } else {
        count = step;
      }
    }
    const bool found = (_processints(it, 0) == idx_configuration);
    return found;
  }

  ordinal_type validateData(const std::vector<std::string> &processes,
                            const std::vector<std::vector<std::vector<ordinal_type>>> &process_constraints) const;
  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;
  std::ostream &showMe(std::ostream &os, const std::string &label, const std::vector<std::string> &processes,
                       const ordinal_type verbose = 0) const;
};

} // namespace KinCat

#endif
