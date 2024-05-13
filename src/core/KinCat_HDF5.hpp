//
//  KinCat_HDF5.h
//  KinCat
//
//  Created by Roberts, Nathan V on 1/4/23, based on EpetraExt_HDF5.
//

#ifndef KinCat_HDF5_h
#define KinCat_HDF5_h

#include "hdf5.h"
#include "KinCat_Config.hpp"

// TODO: initialize HAVE_HDF5 depending on CMake settings.  For now, we hard-code that we have HDF5.
//#define HAVE_HDF5
#ifdef HAVE_HDF5
  #warning "HDF5 Enabled "
#else 
  #warning "HDF5 Not Enabled "
#endif

#include "KinCat_Util.hpp"

#include <string>

namespace KinCat
{

class HDF5
{
public:
  // @{ \name Constructor and destructor.
  //! Constructor
  HDF5();

  //! Destructor
  ~HDF5()
  {
    if (IsOpen())
      Close();
  }

  // @}
  // @{ \name Basic operations
  
  //! Create a new file.
  void Create(const std::string FileName);

  //! Open specified file with given access type.
  void Open(const std::string FileName, int AccessType = H5F_ACC_RDWR);

  //! Close the file.
  void Close()
  {
    H5Fclose(file_id_);
    IsOpen_ = false;
  }
  
  //! Flush the content to the file
  void Flush()
  {
    H5Fflush(file_id_, H5F_SCOPE_GLOBAL);
  }
  
  //! Create group \c GroupName.
  void CreateGroup(const std::string& GroupName)
  {
    hid_t group_id = H5Gcreate(file_id_, GroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
  }
  
  //! Return \c true if \c Name is contained in the database.
  bool IsContained(std::string Name, std::string GroupName = "");
  
  //! Return \c true if a file has already been opened using Open()/Create()
  bool IsOpen() const
  {
    return(IsOpen_);
  }
  
  // @{ \name Distributed arrays
  
  //! Write the distributed array \c data, of type \c type, to group \c GroupName, using dataset name \c DataSetName
  void Write(const std::string& GroupName, const std::string& DataSetName, int MySize, int GlobalSize, hid_t type, const void* data);

  //! Read the distributed array \c data, of type \c type, from group \c GroupName, using dataset name \c DataSetName
  void Read(const std::string& GroupName, const std::string& DataSetName,
            int MySize, int GlobalSize,
            const hid_t type, void* data);

  // @}

  //! Write a double View in group \c GroupName using the given \c DataSetName.
  void WriteDoubleView1D(const std::string& GroupName, const std::string& DataSetName, const value_type_1d_view<double, host_device_type>);
  
  //! Read a double View from group \c GroupName using the given \c DataSetName.
  value_type_1d_view<double, host_device_type> ReadDoubleView1D(const std::string& GroupName, const std::string& DataSetName);
  
  //! Write a char View in group \c GroupName using the given \c DataSetName.
  void WriteCharView1D(const std::string& GroupName, const std::string& DataSetName, const value_type_1d_view<char, host_device_type>);
  
  //! Read a char View from group \c GroupName using the given \c DataSetName.
  value_type_1d_view<char, host_device_type> ReadCharView1D(const std::string& GroupName, const std::string& DataSetName);
  
  //! Write a short View in group \c GroupName using the given \c DataSetName.
  void WriteShortView1D(const std::string& GroupName, const std::string& DataSetName, const value_type_1d_view<short, host_device_type>);
  
  //! Read a short View from group \c GroupName using the given \c DataSetName.
  value_type_1d_view<short, host_device_type> ReadShortView1D(const std::string& GroupName, const std::string& DataSetName);

  //! Write an int View in group \c GroupName using the given \c DataSetName.
  void WriteIntView1D(const std::string& GroupName, const std::string& DataSetName, const value_type_1d_view<int, host_device_type>);
  
  //! Read an int View from group \c GroupName using the given \c DataSetName.
  value_type_1d_view<int, host_device_type> ReadIntView1D(const std::string& GroupName, const std::string& DataSetName);
  
  // @{ \name basic non-distributed data types
  
  //! Write an integer in group \c GroupName using the given \c DataSetName.
  void Write(const std::string& GroupName, const std::string& DataSetName, int data);

  //! Read an integer from group \c /GroupName/DataSetName
  void Read(const std::string& GroupName, const std::string& DataSetName, int& data);

  //! Write a double in group \c GroupName using the given \c DataSetName.
  void Write(const std::string& GroupName, const std::string& DataSetName, double data);

  //! Read a double from group \c /GroupName/DataSetName
  void Read(const std::string& GroupName, const std::string& DataSetName, double& data);

  //! Write a string in group \c GroupName using the given \c DataSetName.
  void Write(const std::string& GroupName, const std::string& DataSetName, const std::string& data);

  //! Read a string from group \c /GroupName/DataSetName
  void Read(const std::string& GroupName, const std::string& DataSetName, std::string& data);

  //! Read the serial array \c data, of type \c type, from group \c GroupName, using the dataset name \c DataSetName.
  void Read(const std::string& GroupName, const std::string& DataSetName,
            const hid_t type, const int Length, void* data);

  //! Write the serial array \c data, of type \c type, to group \c GroupName, using the dataset name \c DataSetName
  void Write(const std::string& GroupName, const std::string& DataSetName,
                       const hid_t type, const int Length,
                       const void* data);

  //! Associate string \c Comment with group \c GroupName.
  void WriteComment(const std::string& GroupName, std::string Comment)
  {
    H5Gset_comment(file_id_, GroupName.c_str(), Comment.c_str());
  }

  //! Read the string associated with group \c GroupName.
  void ReadComment(const std::string& GroupName, std::string& Comment)
  {
    char comment[128];
    H5Gget_comment(file_id_, GroupName.c_str(), 128, comment);
    Comment = comment;
  }

  // @}
  
private:
  // @{ \name Private Data

  //! FileName currently open.
  std::string FileName_;
  //! If \c true, a file is currently open.
  bool IsOpen_;

  //! file ID for HDF5.
  hid_t   file_id_;
  hid_t   plist_id_;
  herr_t  status;

  // @}
};

}

#endif /* KinCat_HDF5_h */
