# Install script for directory: /datalocal/aux0/mnt/patc/Visit_2.10.2_build/src/databases/NETCDF

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libINETCDFDatabase.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libINETCDFDatabase.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libINETCDFDatabase.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases" TYPE SHARED_LIBRARY FILES "/datalocal/aux0/mnt/patc/Visit_2.10.2_build/src/plugins/databases/libINETCDFDatabase.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libINETCDFDatabase.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libINETCDFDatabase.so")
    file(RPATH_REMOVE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libINETCDFDatabase.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libINETCDFDatabase.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libMNETCDFDatabase.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libMNETCDFDatabase.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libMNETCDFDatabase.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases" TYPE SHARED_LIBRARY FILES "/datalocal/aux0/mnt/patc/Visit_2.10.2_build/src/plugins/databases/libMNETCDFDatabase.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libMNETCDFDatabase.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libMNETCDFDatabase.so")
    file(RPATH_REMOVE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libMNETCDFDatabase.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libMNETCDFDatabase.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libENETCDFDatabase_ser.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libENETCDFDatabase_ser.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libENETCDFDatabase_ser.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases" TYPE SHARED_LIBRARY FILES "/datalocal/aux0/mnt/patc/Visit_2.10.2_build/src/plugins/databases/libENETCDFDatabase_ser.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libENETCDFDatabase_ser.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libENETCDFDatabase_ser.so")
    file(RPATH_REMOVE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libENETCDFDatabase_ser.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.10.2/linux-x86_64/plugins/databases/libENETCDFDatabase_ser.so")
    endif()
  endif()
endif()

