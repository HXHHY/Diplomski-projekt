cmake_minimum_required(VERSION 2.8.3)
project(ceres_mpc)

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

find_package(catkin REQUIRED COMPONENTS
    roscpp)

# reading usr/lib/cmake/eigen3/Eigen3Config.cmake
find_package(Eigen3 REQUIRED)
# reading usr/local/lib/cmake/Ceres/CeresConfig.cmake
find_package(Ceres REQUIRED)

# DEBUG output
#MESSAGE(STATUS "catkin include dirs: 
#    ${catkin_INCLUDE_DIRS}")
#MESSAGE(STATUS "eigen include dirs: 
#    ${EIGEN3_INCLUDE_DIRS}")

catkin_package(
    INCLUDE_DIRS
        include
        ${catkin_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIRS}
    CATKIN_DEPENDS
        roscpp
    DEPENDS
        cmake_modules
        Eigen3
        Ceres
    LIBRARIES
        mpc_mm_ceres
        ${catkin_LIBRARIES}
        ${CERES_LIBRARIES})

include_directories(
   include
   ${catkin_INCLUDE_DIRS}
   ${EIGEN3_INCLUDE_DIRS}
)

# EIGEN_LIBRARIES NE POSTOJE !
add_executable(ceres_mpc_1 src/ceres_mpc_1.cpp)
target_link_libraries(ceres_mpc_1 
     ${catkin_LIBRARIES} 
     ${CERES_LIBRARIES})

add_executable(ceres_mpc_5 src/ceres_mpc_5.cpp)
target_link_libraries(ceres_mpc_5 
    ${catkin_LIBRARIES} 
    ${CERES_LIBRARIES})

add_executable(ceres_mpc_4 src/ceres_mpc_4.cpp)
target_link_libraries(ceres_mpc_4 
    ${catkin_LIBRARIES} 
    ${CERES_LIBRARIES})

# library
add_library(mpc_mm_ceres
        src/mpc_mm_ceres.cpp)
target_link_libraries(mpc_mm_ceres
    ${catkin_LIBRARIES} 
    ${CERES_LIBRARIES})

# Executable za probu
add_executable(ceres_proba_razreda 
    src/proba_mpc_ceres_razreda.cpp)
target_link_libraries(ceres_proba_razreda mpc_mm_ceres)
