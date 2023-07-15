cmake_minimum_required(VERSION 3.0.2)
project(test_ros)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
)
find_package(PCL REQUIRED)
find_package( Threads )

file(GLOB_RECURSE SRC_DIR_LIST "src/*.cpp")	
include_directories(
# include
  ${catkin_INCLUDE_DIRS}
)
include_directories(${PCL_INCLUDE_DIRS})

aux_source_directory(./src SRCS)
add_executable(test_ros ${SRC_DIR_LIST}	)
target_link_libraries(test_ros ${catkin_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${PCL_LIBRARIES})
