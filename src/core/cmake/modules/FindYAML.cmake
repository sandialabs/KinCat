# Check yaml installation
FILE(GLOB_RECURSE YAML_FOUND_CMAKE_FILE "${YAML_INSTALL_PATH}/yaml-cpp-config.cmake")

IF (YAML_FOUND_CMAKE_FILE)
  
  INCLUDE(${YAML_FOUND_CMAKE_FILE})

  # when cmake includes YAML_CPP_INCLUDE_DIR variable
  IF (YAML_CPP_INCLUDE_DIR)
    SET(KINCAT_YAML_INCLUDE_PATH ${YAML_CPP_INCLUDE_DIR})
  ELSE()
    GET_TARGET_PROPERTY(KINCAT_YAML_INCLUDE_PATH yaml-cpp INTERFACE_INCLUDE_DIRECTORIES)
  ENDIF()

  # when cmake includes YAML_CPP_LIBRARIES
  GET_TARGET_PROPERTY(KINCAT_YAML_LIBRARIES_IMPORTED_CONFIGURATION yaml-cpp IMPORTED_CONFIGURATIONS)
  SET(KINCAT_YAML_LIBRARIES_PROPERTY_NAME "IMPORTED_LOCATION_${KINCAT_YAML_LIBRARIES_IMPORTED_CONFIGURATION}")
  GET_TARGET_PROPERTY(KINCAT_YAML_LIBRARIES yaml-cpp ${KINCAT_YAML_LIBRARIES_PROPERTY_NAME})
  
  MESSAGE(STATUS "  YAML include dir ${KINCAT_YAML_INCLUDE_PATH}")
  MESSAGE(STATUS "  YAML library ${KINCAT_YAML_LIBRARIES}")  
  SET(YAML_FOUND ON)
  MESSAGE("-- YAML is found at ${YAML_INSTALL_PATH}")  
  
  ADD_LIBRARY(yaml INTERFACE)
  SET_TARGET_PROPERTIES(yaml PROPERTIES 
    INTERFACE_INCLUDE_DIRECTORIES "${KINCAT_YAML_INCLUDE_PATH}"
    INTERFACE_COMPILE_OPTIONS "-I${KINCAT_YAML_INCLUDE_PATH}"
    INTERFACE_LINK_LIBRARIES "${KINCAT_YAML_LIBRARIES}"
  )    
ELSE()
  MESSAGE(FATAL_ERROR "-- YAML is not found at ${YAML_INSTALL_PATH}")
ENDIF()
