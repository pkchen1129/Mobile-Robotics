cmake_minimum_required(VERSION 3.9.1)
project(Mobile_ps4)

##### Compile as C++11 #####
add_compile_options(-std=c++11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)  #set the output path to /bin


#### System dependencies are found with CMake's conventions ####
find_package(GTSAM REQUIRED)
find_package(Eigen3 REQUIRED)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/dataset
  ${GTSAM_INCLUDE_DIR}
  ${EIGEN3_INCLUDE_DIR}
)
# message(STATUS ${GTSAM_INCLUDE_DIR})

## Declare a C++ library
add_library( lib ${PROJECT_SOURCE_DIR}/src/main.cpp )
add_executable(test src/main.cpp)
target_link_libraries( 
   test
   lib
   gtsam
   ${EIGEN3_LIBS}
   )
