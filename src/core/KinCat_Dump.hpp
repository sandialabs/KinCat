#ifndef __KINCAT_DUMP_HPP__
#define __KINCAT_DUMP_HPP__

#include "KinCat_Lattice.hpp"
#include "KinCat_Util.hpp"

namespace KinCat {

template <typename DeviceType> struct Dump {
private:
  std::string _filename_dump, _filename_site;
  Lattice<DeviceType> _lattice;

  std::ofstream _ofs;
  Lattice<host_device_type> _lattice_host;

  bool _is_first;

  ordinal_type snapshot(const ordinal_type sid, const real_type t,
                        const value_type_1d_view<site_type, host_device_type> sites, const ordinal_type verbose = 0);

public:
  Dump() = delete;
  Dump(const Dump &b) = delete;
  Dump(const std::string &filename_dump, const std::string &filename_site, const Lattice<DeviceType> &lattice);
  virtual ~Dump() = default;

  ordinal_type initialize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);
  ordinal_type finalize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);
  ordinal_type snapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);

  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;
};

} // namespace KinCat

#endif
