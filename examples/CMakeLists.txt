# 创建示例可执行文件
add_executable(ckks_demo ckks_demo.cpp)

target_include_directories(ckks_demo PRIVATE 
    ${CMAKE_SOURCE_DIR}/include 
    ${OPENFHE_INCLUDE_DIRS}
    )
target_link_libraries(ckks_demo ${OpenFHE_LIBRARIES})


add_executable(opfhe_trapdoor_demo openfhe_trapdoortest.cpp)
target_include_directories(opfhe_trapdoor_demo PRIVATE 
    ${CMAKE_SOURCE_DIR}/include 
    ${OPENFHE_INCLUDE_DIRS}
    )
target_link_libraries(opfhe_trapdoor_demo ${OpenFHE_LIBRARIES})