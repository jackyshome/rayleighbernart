cmake_minimum_required (VERSION 2.8.1)

project (rayleighbernart)

## #############################################################################
## cmake module path
## #############################################################################

set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/find
    ${CMAKE_MODULE_PATH}
    )

## #############################################################################
## enable c++ 11
## #############################################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## #############################################################################
## find third parts
## #############################################################################


find_package(CImg REQUIRED)
find_package(FFTW REQUIRED)
find_package(smctc REQUIRED)
find_package(X11 REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)


## #############################################################################
## include directories
## #############################################################################

include_directories(include
    ${CImg_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
    ${smctc_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    /opt/X11/include
    )

link_directories(${smctc_LIBRARY_DIRS}
    /opt/X11/lib
    )

## #############################################################################
## source defenition
## #############################################################################

set(${PROJECT_NAME}_CPP
  DynamicalModel.cpp  
  main.cpp  
  ObservationModel.cpp  
  PffuncsRayleigh.cpp  
  ReducedOrderModel.cpp  
  )

set(${PROJECT_NAME}_H 
  include/DynamicalModel.h
  include/ObservationModel.h
  include/PffuncsRayleigh.h
  include/ReducedOrderModel.h
  )

## #############################################################################
## add executable
## #############################################################################

add_executable(rayleighbernart
  ${${PROJECT_NAME}_CPP}
  ${${PROJECT_NAME}_H}
  )

add_executable(DNS
  DNS.cpp
  DynamicalModel.cpp
  ObservationModel.cpp
  )

## #############################################################################
## link third parts
## #############################################################################

message("smctc_LIBRARY_DIRS: ${smctc_LIBRARY_DIRS}")

target_link_libraries(${PROJECT_NAME}
    smctc
    X11
    fftw3
    pthread
    gsl
    lapack
    boost_system
    boost_filesystem
    )

target_link_libraries(DNS
    X11
    pthread
    fftw3
    gsl
    lapack
    )
