#ifndef __KINCAT_DUMP_HPP__
#define __KINCAT_DUMP_HPP__

#include "KinCat_Lattice.hpp"
#include "KinCat_Util.hpp"

#include "KinCat_HDF5.hpp"

namespace KinCat {

template <typename DeviceType> struct Dump {
private:
  std::string _filename_dump, _filename_site;
  Lattice<DeviceType> _lattice;

  std::ofstream _ofs;
  Lattice<host_device_type> _lattice_host;

  bool _is_first;

  ordinal_type _nextSnapshotID = 0;
  
  ordinal_type snapshot(const ordinal_type sid, const real_type t,
                        const value_type_1d_view<site_type, host_device_type> sites, const bool restart_flag, const ordinal_type verbose = 0);

#ifdef HAVE_HDF5
  HDF5 _hdf5;
#endif
  bool _useHDF5;
public:
  Dump() = delete;
  Dump(const Dump &b) = delete;
  Dump(const std::string &filename_dump, const std::string &filename_site, const Lattice<DeviceType> &lattice,
       const bool useHDF5 = false);
  virtual ~Dump() = default;

  ordinal_type initialize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);
  ordinal_type finalize(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);
  ordinal_type snapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);
  ordinal_type sitessnapshot(const value_type_1d_view<real_type, host_device_type> t, const ordinal_type verbose = 0);

  virtual std::ostream &showMe(std::ostream &os, const std::string &label, const ordinal_type verbose = 0) const;
};

} // namespace KinCat

#endif