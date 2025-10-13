#include <iostream>
#include "util.h"
#include <openfhe.h>
// main.cpp
#include "MyRLWETrapdoor.h"

using namespace lbcrypto;
using namespace myrlwe;

template <class Element>
void demonstrateMyRLWE() {
    // 1. 初始化OpenFHE上下文
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 50;
    uint32_t batchSize = 8;
    SecurityLevel securityLevel = HEStd_128_classic;
    
    // 创建加密参数
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleFactorBits);
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(securityLevel);

    
    CryptoContext<Element> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    
    // 2. 生成密钥
    // auto keys = cc->KeyGen();
    // cc->EvalMultKeyGen(keys.secretKey);
    
    // 3. 获取环参数
    auto params = cc->GetElementParams();
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    auto uniform_alloc = Element::MakeDiscreteUniformAllocator(params, EVALUATION);
        std::cout << "开始在外面生成g"<< std::endl;
    Matrix<Element> g = Matrix<Element>(zero_alloc, 1, 5).GadgetVector(2);
    // 4. 创建a'矩阵 (示例：1×2矩阵)
    Matrix<Element> a_prime(zero_alloc, 1, 64);
    a_prime(0, 0) = uniform_alloc();
    a_prime(0, 1) = uniform_alloc();
    
    // 5. 创建h元素
    Element h = uniform_alloc();

    // 6. 使用我们的Trapdoor生成算法
    double stddev = 4.0;
    auto trapdoor_result = MyRLWETrapdoorUtility<Element>::TrapdoorGen(
        params, a_prime, h, stddev);
    
    Matrix<Element> A = trapdoor_result.first;
    MyRLWETrapdoorPair<Element> T_pair = trapdoor_result.second;
    
    std::cout << "Trapdoor generation completed!" << std::endl;
    
    // 7. 使用SamplePre算法
    Element u = uniform_alloc();  // 目标syndrome
    double eta = 8.0, sigma = 3.2, alpha = 1.0;
    
    typename Element::DggType dgg(stddev);
    
    Matrix<Element> x = MyRLWETrapdoorUtility<Element>::SamplePre(
        T_pair, A, h, eta, sigma, alpha, u, dgg);
    
    std::cout << "SamplePre completed!" << std::endl;
    
    // 8. 验证: 检查 A·x ≈ h·u
    Matrix<Element> Ax = A * x;
    Element hu = h * u;
    
    // 计算误差
    Element error = Ax(0, 0) - hu;
    error.SetFormat(COEFFICIENT);
    
    std::cout << "Verification completed!" << std::endl;
    // 可以进一步分析error的大小来验证正确性
}

int main() {
    std::cout << "trapdooribe fjj" << std::endl;
    std::cout << "================" << std::endl;
    std::cout << "Testing with Poly..." << std::endl;
    demonstrateMyRLWE<Poly>();
    // std::cout << "Testing with DCRTPoly..." << std::endl;
    // demonstrateMyRLWE<DCRTPoly>();
    return 0;
}
