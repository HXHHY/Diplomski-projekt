cmake_minimum_required(VERSION 2.8.3)
project(ceres_mpc)



if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

find_package(catkin REQUIRED
            COMPONENTS
            roscpp)
find_package(cmake_modules REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)

include_directories(${catkin_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CERES_INCLUDE_DIRS})

add_definitions(${EIGEN_DEFINITIONS})

catkin_package(
    INCLUDE_DIRS
        include
        ${catkin_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${CERES_INCLUDE_DIRS}
    CATKIN_DEPENDS
        roscpp
    DEPENDS
        cmake_modules
        Eigen3
        Ceres)
include_directories(
   include
   ${catkin_INCLUDE_DIRS}
)
add_executable(ceres_mpc_1 src/ceres_mpc_1.cpp)

target_link_libraries(ceres_mpc_1 ${catkin_LIBRARIES} ${CERES_LIBRARIES} ${EIGEN_LIBRARIES})

add_executable(ceres_mpc_5 src/ceres_mpc_5.cpp)

target_link_libraries(ceres_mpc_5 ${catkin_LIBRARIES} ${CERES_LIBRARIES} ${EIGEN_LIBRARIES})

add_executable(ceres_mpc_4 src/ceres_mpc_4.cpp)

add_executable(ceres_proba_razreda src/proba_mpc_ceres_razreda.cpp)

target_link_libraries(ceres_mpc_4 ${catkin_LIBRARIES} ${CERES_LIBRARIES} ${EIGEN_LIBRARIES})

add_library(mpc_mm_ceres
        src/mpc_mm_ceres.cpp
        include/ceres_mpc/mpc_mm_ceres.h

)

target_link_libraries(mpc_mm_ceres
 ${catkin_LIBRARIES} ${CERES_LIBRARIES} ${EIGEN_LIBRARIES}
)

target_link_libraries(ceres_proba_razreda mpc_mm_ceres)
