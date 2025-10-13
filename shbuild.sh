rm -rf build
mkdir build
cd build

cmake ..

make -j2

# # 运行CKKS示例
# ./examples/ckks_demo
# ./examples/opfhe_trapdoor_demo

./src/trapibe_openfhe_fjj
