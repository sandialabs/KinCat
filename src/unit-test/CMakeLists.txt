#
# Craig test
#
ADD_EXECUTABLE(kincat-test.x KinCat_Test_Driver.cpp)
TARGET_LINK_LIBRARIES(kincat-test.x ${KINCAT_LINK_LIBRARIES})

INSTALL(TARGETS kincat-test.x
        PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE
        DESTINATION "${CMAKE_INSTALL_PREFIX}/${KINCAT_INSTALL_TEST_PATH}")

ADD_TEST(KinCat_UnitTest kincat-test.x)

#
# Testing files
#
FILE(COPY test-files
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

INSTALL(DIRECTORY test-files
        DESTINATION ${CMAKE_INSTALL_PREFIX}/unit-test/test-files)
