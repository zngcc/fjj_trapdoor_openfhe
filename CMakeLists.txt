cmake_minimum_required(VERSION 3.15)
project(OpenFHEDemo VERSION 1.0.0 LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找OpenFHE库
set(OpenFHE_DIR "/usr/local/lib/OpenFHE")
set(OPENFHE_INCLUDE_DIRS
    "/usr/local/include/openfhe"
    "/usr/local/include/openfhe/core"
    "/usr/local/include/openfhe/pke"
    "/usr/local/include/openfhe/binfhe"
    "/usr/local/include/openfhe/cereal"
)
find_package(OpenFHE REQUIRED)
# 包含目录
include_directories(include)
include_directories(${OpenFHE_INCLUDE_DIRS})
# 添加子目录
add_subdirectory(src)

# 如果存在examples目录，也添加
if(EXISTS "${CMAKE_SOURCE_DIR}/examples")
    add_subdirectory(examples)
endif()