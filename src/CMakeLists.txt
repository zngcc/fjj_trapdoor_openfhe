# 创建主可执行文件
add_executable(trapibe_openfhe_fjj
    main.cpp
    MyRLWETrapdoor.cpp
)

# 设置目标属性
target_include_directories(trapibe_openfhe_fjj PRIVATE 
    ${CMAKE_SOURCE_DIR}/include 
    ${OPENFHE_INCLUDE_DIRS}
    )

# 链接OpenFHE库
target_link_libraries(trapibe_openfhe_fjj ${OpenFHE_LIBRARIES})
target_compile_features(trapibe_openfhe_fjj PRIVATE cxx_std_17)