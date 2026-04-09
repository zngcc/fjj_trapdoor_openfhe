# 创建示例可执行文件
add_executable(ckks_demo ckks_demo.cpp)

target_include_directories(ckks_demo PRIVATE 
    ${CMAKE_SOURCE_DIR}/include 
    ${OPENFHE_INCLUDE_DIRS}
    )
target_link_libraries(ckks_demo ${OpenFHE_LIBRARIES})

add_executable(dcrt_demo dcrt_demo.cpp)
target_include_directories(dcrt_demo PRIVATE 
    ${CMAKE_SOURCE_DIR}/include 
    ${OPENFHE_INCLUDE_DIRS}
    )
target_link_libraries(dcrt_demo ${OpenFHE_LIBRARIES})

