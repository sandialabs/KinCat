//
//  KinCat_HDF5.cpp
//  KinCat
//
//  Created by Roberts, Nathan V on 1/4/23.
//
// Most code copied from EpetraExt_HDF5.cpp.

#include "KinCat_HDF5.hpp"
#include "KinCat_Util.hpp"

using namespace KinCat;

namespace {
// ==========================================================================
// data container and iterators to find a dataset with a given name
struct FindDataset_t {
  std::string name;
  bool found;
};

static herr_t FindDataset(hid_t loc_id, const char *name, void *opdata) {
  std::string &token = ((FindDataset_t *)opdata)->name;
  if (token == name)
    ((FindDataset_t *)opdata)->found = true;

  return (0);
}
} // namespace

HDF5::HDF5() : IsOpen_(false) {}

void HDF5::Create(const std::string FileName) {
  KINCAT_CHECK_ERROR(
      IsOpen(),
      "an HDF5 is already open; first close the current one using method Close(), then open/create a new one");

  FileName_ = FileName;

  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);

  // create the file collectively and release property list identifier.
  file_id_ = H5Fcreate(FileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id_);
  H5Pclose(plist_id_);

  IsOpen_ = true;
}

// ==========================================================================
bool HDF5::IsContained(std::string Name, std::string GroupName) {
  KINCAT_CHECK_ERROR(!IsOpen(), "no file open yet");

  FindDataset_t data;
  data.name = Name;
  data.found = false;

  // recursively look for groups
  size_t pos = Name.find("/");
  if (pos != std::string::npos) {
    std::string NewGroupName = Name.substr(0, pos);
    if (GroupName != "")
      NewGroupName = GroupName + "/" + NewGroupName;
    std::string NewName = Name.substr(pos + 1);
    return IsContained(NewName, NewGroupName);
  }

  GroupName = "/" + GroupName;

  // int idx_f =
  H5Giterate(file_id_, GroupName.c_str(), NULL, FindDataset, (void *)&data);

  return (data.found);
}

void HDF5::Open(const std::string FileName, int AccessType) {
  KINCAT_CHECK_ERROR(
      IsOpen(), "an HDF5 is already open; first close the current one using method Close(), then open/create a new one")

  FileName_ = FileName;

  // Set up file access property list with parallel I/O access
  plist_id_ = H5Pcreate(H5P_FILE_ACCESS);

  // create the file collectively and release property list identifier.
  file_id_ = H5Fopen(FileName.c_str(), AccessType, plist_id_);
  H5Pclose(plist_id_);

  IsOpen_ = true;
}

// ==========================================================================
void HDF5::Write(const std::string &GroupName, const std::string &DataSetName, int what) {
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id =
      H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_INT, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id, H5P_DEFAULT, &what);
  KINCAT_CHECK_ERROR((status < 0), "function H5Dwrite returned a negative value");

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void HDF5::Write(const std::string &GroupName, const std::string &DataSetName, double what) {
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id =
      H5Dcreate(group_id, DataSetName.c_str(), H5T_NATIVE_DOUBLE, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_id, H5P_DEFAULT, &what);
  KINCAT_CHECK_ERROR((status < 0), "function H5Giterater returned a negative value");

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void HDF5::Read(const std::string &GroupName, const std::string &DataSetName, int &data) {
  KINCAT_CHECK_ERROR(!IsContained(GroupName), "requested group not found");

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, filespace_id, H5P_DEFAULT, &data);
  KINCAT_CHECK_ERROR((status < 0), "function H5Giterater returned a negative value");

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

// ==========================================================================
void HDF5::Read(const std::string &GroupName, const std::string &DataSetName, double &data) {
  KINCAT_CHECK_ERROR(!IsContained(GroupName), "requested group not found");

  // Create group in the root group using absolute name.
  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t filespace_id = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace_id, H5P_DEFAULT, &data);
  KINCAT_CHECK_ERROR((status < 0), "function H5Giterater returned a negative value");

  H5Sclose(filespace_id);
  H5Dclose(dset_id);
  H5Gclose(group_id);
}

// ==========================================================================
void HDF5::Write(const std::string &GroupName, const std::string &DataSetName, const std::string &data) {
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hsize_t len = 1;

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataspace_id = H5Screate_simple(1, &len, NULL);

  hid_t atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, data.size() + 1);

  hid_t dataset_id =
      H5Dcreate(group_id, DataSetName.c_str(), atype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dataset_id, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.c_str());
  KINCAT_CHECK_ERROR((status < 0), "function H5Dwrite returned a negative value");

  status = H5Tclose(atype);
  KINCAT_CHECK_ERROR((status < 0), "function H5Tclose returned a negative value");
  status = H5Dclose(dataset_id);
  KINCAT_CHECK_ERROR((status < 0), "function H5Dclose returned a negative value");
  status = H5Sclose(dataspace_id);
  KINCAT_CHECK_ERROR((status < 0), "function H5Sclose returned a negative value");
  status = H5Gclose(group_id);
  KINCAT_CHECK_ERROR((status < 0), "function H5Gclose returned a negative value");
}

// ==========================================================================
void HDF5::Read(const std::string &GroupName, const std::string &DataSetName, std::string &data) {
  KINCAT_CHECK_ERROR(!IsContained(GroupName), "requested group not found");

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  hid_t datatype_id = H5Dget_type(dataset_id);
  //  size_t typesize_id = H5Tget_size(datatype_id);
  H5T_class_t typeclass_id = H5Tget_class(datatype_id);

  KINCAT_CHECK_ERROR(typeclass_id != H5T_STRING, "requested group is not a std::string");

  char data2[160];
  herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2);
  KINCAT_CHECK_ERROR((status < 0), "function H5Dread returned a negative value");

  data = data2;

  H5Dclose(dataset_id);
  H5Gclose(group_id);
}

// ============= //
// serial arrays //
// ============= //

// ==========================================================================
void HDF5::Write(const std::string &GroupName, const std::string &DataSetName, const hid_t type, const int Length,
                 const void *data) {
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hsize_t dimsf = Length;

  hid_t filespace_id = H5Screate_simple(1, &dimsf, NULL);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), type, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t status = H5Dwrite(dset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  KINCAT_CHECK_ERROR((status < 0), "function H5Dwrite returned a negative value");

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ==========================================================================
void HDF5::Read(const std::string &GroupName, const std::string &DataSetName, const hid_t type, const int Length,
                void *data) {
  KINCAT_CHECK_ERROR(!IsContained(GroupName), "requested group not found")

  hsize_t dimsf = Length;

  hid_t filespace_id = H5Screate_simple(1, &dimsf, NULL);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);

  hid_t dset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);

  herr_t status = H5Dread(dset_id, type, H5S_ALL, filespace_id, H5P_DEFAULT, data);
  KINCAT_CHECK_ERROR((status < 0), "function H5Dread returned a negative value");

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
}

// ================== //
// distributed arrays //
// ================== //

// ==========================================================================
void HDF5::Write(const std::string &GroupName, const std::string &DataSetName, int MySize, int GlobalSize, hid_t type,
                 const void *data) {
  int Offset = 0;

  hsize_t MySize_t = MySize;
  hsize_t GlobalSize_t = GlobalSize;
  hsize_t Offset_t = Offset;

  hid_t filespace_id = H5Screate_simple(1, &GlobalSize_t, NULL);

  // Create the dataset with default properties and close filespace.
  if (!IsContained(GroupName))
    CreateGroup(GroupName);

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dset_id = H5Dcreate(group_id, DataSetName.c_str(), type, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Sclose(filespace_id);

  // Each process defines dataset in memory and writes it to the hyperslab
  // in the file.

  hid_t memspace_id = H5Screate_simple(1, &MySize_t, NULL);

  // Select hyperslab in the file.
  filespace_id = H5Dget_space(dset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &Offset_t, NULL, &MySize_t, NULL);

  plist_id_ = H5Pcreate(H5P_DATASET_XFER);

  status = H5Dwrite(dset_id, type, memspace_id, filespace_id, plist_id_, data);
  KINCAT_CHECK_ERROR((status < 0), "function H5Dwrite returned a negative value");

  // Close/release resources.
  H5Dclose(dset_id);
  H5Gclose(group_id);
  H5Sclose(filespace_id);
  H5Sclose(memspace_id);
  H5Pclose(plist_id_);
}

// ==========================================================================
void HDF5::Read(const std::string &GroupName, const std::string &DataSetName, int MySize, int GlobalSize,
                const hid_t type, void *data) {
  KINCAT_CHECK_ERROR(!IsOpen(), "no file open yet");

  hsize_t MySize_t = MySize;

  // offset
  hsize_t Offset_t = 0;

  hid_t group_id = H5Gopen(file_id_, GroupName.c_str(), H5P_DEFAULT);
  hid_t dataset_id = H5Dopen(group_id, DataSetName.c_str(), H5P_DEFAULT);
  // hid_t space_id = H5Screate_simple(1, &Offset_t, 0);

  // Select hyperslab in the file.
  hid_t filespace_id = H5Dget_space(dataset_id);
  H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, &Offset_t, NULL, &MySize_t, NULL);

  hid_t mem_dataspace = H5Screate_simple(1, &MySize_t, NULL);

  herr_t status = H5Dread(dataset_id, type, mem_dataspace, filespace_id, H5P_DEFAULT, data);
  KINCAT_CHECK_ERROR((status < 0), "function H5Dread returned a negative value");

  H5Sclose(mem_dataspace);
  H5Gclose(group_id);
  // H5Sclose(space_id);
  H5Dclose(dataset_id);
  //  H5Dclose(filespace_id);
}

value_type_1d_view<char, host_device_type> HDF5::ReadCharView1D(const std::string &GroupName, const std::string &DataSetName) {
  // read the length first
  std::string DataSetLengthName = DataSetName + ".Length";
  
  int length;
  Read(GroupName, DataSetLengthName, length);
  
  value_type_1d_view<char, host_device_type> dataView(GroupName + "/" + DataSetName, length);
  
  hid_t type = H5T_NATIVE_CHAR;
  char *view_raw = dataView.data(); // Discouraged practice. Should be fine for 1D views
  Read(GroupName, DataSetName, length, type, view_raw);
  
  return dataView;
}

void HDF5::WriteShortView1D(const std::string &GroupName, const std::string &DataSetName,
                           const value_type_1d_view<short, host_device_type> dataView) {
  hid_t type = H5T_NATIVE_SHORT;

  const short *view_raw = dataView.data(); // Discouraged practice. Should be fine for 1D views
  const int length = dataView.extent(0);
  
  std::string DataSetLengthName = DataSetName + ".Length";
  Write(GroupName, DataSetLengthName, length);
  Write(GroupName, DataSetName, type, length, view_raw);
}

void HDF5::WriteCharView1D(const std::string &GroupName, const std::string &DataSetName,
                           const value_type_1d_view<char, host_device_type> dataView) {
  hid_t type = H5T_NATIVE_CHAR;

  const char *view_raw = dataView.data(); // Discouraged practice. Should be fine for 1D views
  const int length = dataView.extent(0);
  
  std::string DataSetLengthName = DataSetName + ".Length";
  Write(GroupName, DataSetLengthName, length);
  Write(GroupName, DataSetName, type, length, view_raw);
}

value_type_1d_view<int, host_device_type> HDF5::ReadIntView1D(const std::string &GroupName, const std::string &DataSetName) {
  // read the length first
  std::string DataSetLengthName = DataSetName + ".Length";
  
  int length;
  Read(GroupName, DataSetLengthName, length);
  
  value_type_1d_view<int, host_device_type> dataView(GroupName + "/" + DataSetName, length);
  
  hid_t type = H5T_NATIVE_INT;
  int *view_raw = dataView.data(); // Discouraged practice. Should be fine for 1D views
  Read(GroupName, DataSetName, length, type, view_raw);
  
  return dataView;
}

void HDF5::WriteIntView1D(const std::string &GroupName, const std::string &DataSetName,
                          const value_type_1d_view<int, host_device_type> dataView) {
  hid_t type = H5T_NATIVE_INT;

  const int *view_raw = dataView.data(); // Discouraged practice. Should be fine for 1D views
  const int length = dataView.extent(0);
  
  std::string DataSetLengthName = DataSetName + ".Length";
  Write(GroupName, DataSetLengthName, length);
  Write(GroupName, DataSetName, type, length, view_raw);
}

value_type_1d_view<double, host_device_type> HDF5::ReadDoubleView1D(const std::string &GroupName, const std::string &DataSetName) {
  // read the length first
  std::string DataSetLengthName = DataSetName + ".Length";
  
  int length;
  Read(GroupName, DataSetLengthName, length);
  
  value_type_1d_view<double, host_device_type> dataView(GroupName + "/" + DataSetName, length);
  
  hid_t type = H5T_NATIVE_DOUBLE;
  double *view_raw = dataView.data(); // Discouraged practice. Should be fine for 1D views
  Read(GroupName, DataSetName, length, type, view_raw);
  
  return dataView;
}

void HDF5::WriteDoubleView1D(const std::string &GroupName, const std::string &DataSetName,
                             const value_type_1d_view<double, host_device_type> dataView) {
  hid_t type = H5T_NATIVE_DOUBLE;

  const double *view_raw = dataView.data(); // Discouraged practice. Should be fine for 1D views
  const int length = dataView.extent(0);
  
  std::string DataSetLengthName = DataSetName + ".Length";
  Write(GroupName, DataSetLengthName, length);
  Write(GroupName, DataSetName, type, length, view_raw);
}
