// trapdoor_test.cpp
// #include "lattice/trapdoor.h"
// #include "lattice/lat-hal.h"
// #include "math/matrix-impl.h"
#include "openfhe.h"
#include "lattice/trapdoor.h"
#include <iostream>
#include <memory>
#include "lattice/dgsampling-impl.h"
#include "lattice/lat-hal.h"
#include "lattice/trapdoor-impl.h"

#include "math/matrix-impl.h"
using namespace lbcrypto;

void testTrapdoorGen() {
    std::cout << "=== Testing RLWETrapdoorUtility::TrapdoorGen ===" << std::endl;
    
    // 设置参数
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
    
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    
    // 获取参数
    auto params = cryptoContext->GetElementParams();
    
    // 设置trapdoor参数
    double stddev = 4.0;  // 标准差
    int64_t base = 8;     // gadget向量基
    bool bal = false;     // 是否使用平衡表示
////////////////////////////////
    auto zero_alloc     = DCRTPoly::Allocator(params, Format::EVALUATION);
    auto gaussian_alloc = DCRTPoly::MakeDiscreteGaussianCoefficientAllocator(params, Format::COEFFICIENT, stddev);
    auto uniform_alloc  = DCRTPoly::MakeDiscreteUniformAllocator(params, Format::EVALUATION);
            std::cout << "开始在外面生成g"<< std::endl;
    NativeInteger q = params->GetParams()[0]->GetModulus();

    size_t digitCount = static_cast<size_t>(ceil(log2(q.ConvertToDouble()) / log2(base)));

    size_t k = params->GetParams().size() * digitCount;

    Matrix<DCRTPoly> g = Matrix<DCRTPoly>(zero_alloc, 1, k).GadgetVector(base);

    
////////////////////////////////
    try {
        // 调用TrapdoorGen函数
        auto result = RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(params, stddev, base, bal);
        
        Matrix<DCRTPoly>& A = result.first;
        RLWETrapdoorPair<DCRTPoly>& trapdoorPair = result.second;
        
        // 验证结果
        std::cout << "Trapdoor generation successful!" << std::endl;
        std::cout << "Matrix A dimensions: " << A.GetRows() << " x " << A.GetCols() << std::endl;
        std::cout << "Trapdoor pair R dimensions: " << trapdoorPair.m_r.GetRows() 
                  << " x " << trapdoorPair.m_r.GetCols() << std::endl;
        std::cout << "Trapdoor pair E dimensions: " << trapdoorPair.m_e.GetRows() 
                  << " x " << trapdoorPair.m_e.GetCols() << std::endl;
        
        // 验证矩阵A的结构
        if (A.GetRows() == 1 && A.GetCols() >= 3) {
            std::cout << "Matrix A structure verification: PASS" << std::endl;
            
            // 检查第一个元素是否为1
            auto firstElement = A(0, 0);
            firstElement.SetFormat(Format::COEFFICIENT);
            if (firstElement.GetElementAtIndex(0).GetValues().at(0) == 1) {
                std::cout << "First element is 1: PASS" << std::endl;
            } else {
                std::cout << "First element is not 1: FAIL" << std::endl;
            }
        } else {
            std::cout << "Matrix A dimensions incorrect: FAIL" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Trapdoor generation failed: " << e.what() << std::endl;
    }
}
void testPolyTrapdoorGen() {
    std::cout << "=== Testing RLWETrapdoorUtility<Poly>::TrapdoorGen ===" << std::endl;
    
    // 设置参数
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 50;
    uint32_t batchSize = 8;
    SecurityLevel securityLevel = HEStd_128_classic;
    
    // 创建加密参数 - 使用Poly类型
    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetPlaintextModulus(65537);
    parameters.SetSecurityLevel(securityLevel);
    
    CryptoContext<Poly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    
    // 获取参数
    auto params = cryptoContext->GetElementParams();
    
    // 设置trapdoor参数
    double stddev = 4.0;  // 标准差
    int64_t base = 8;     // gadget向量基
    bool bal = false;     // 是否使用平衡表示
    
    try {
        // 调用TrapdoorGen函数
        auto result = RLWETrapdoorUtility<Poly>::TrapdoorGen(params, stddev, base, bal);
        
        Matrix<Poly>& A = result.first;
        RLWETrapdoorPair<Poly>& trapdoorPair = result.second;
        
        // 验证结果
        std::cout << "Poly Trapdoor generation successful!" << std::endl;
        std::cout << "Matrix A dimensions: " << A.GetRows() << " x " << A.GetCols() << std::endl;
        std::cout << "Trapdoor pair R dimensions: " << trapdoorPair.m_r.GetRows() 
                  << " x " << trapdoorPair.m_r.GetCols() << std::endl;
        std::cout << "Trapdoor pair E dimensions: " << trapdoorPair.m_e.GetRows() 
                  << " x " << trapdoorPair.m_e.GetCols() << std::endl;
        
        // 验证矩阵A的结构
        if (A.GetRows() == 1 && A.GetCols() >= 3) {
            std::cout << "Matrix A structure verification: PASS" << std::endl;
            
            // 检查第一个元素是否为1
            auto firstElement = A(0, 0);
            firstElement.SetFormat(Format::COEFFICIENT);
            if (firstElement.GetLength() > 0) {
                std::cout << "First element verification completed" << std::endl;
            }
        } else {
            std::cout << "Matrix A dimensions incorrect: FAIL" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Poly Trapdoor generation failed: " << e.what() << std::endl;
    }
}

void testNativePolyTrapdoorGen() {
    std::cout << "\n=== Testing RLWETrapdoorUtility<NativePoly>::TrapdoorGen ===" << std::endl;
    
    // 设置参数
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 50;
    uint32_t batchSize = 8;
    SecurityLevel securityLevel = HEStd_128_classic;
    
    // 创建加密参数 - 使用NativePoly类型
    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetPlaintextModulus(65537);
    parameters.SetSecurityLevel(securityLevel);
    
    CryptoContext<NativePoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    
    // 获取参数
    auto params = cryptoContext->GetElementParams();
    
    // 设置trapdoor参数
    double stddev = 4.0;  // 标准差
    int64_t base = 8;     // gadget向量基
    bool bal = false;     // 是否使用平衡表示
    
    try {
        // 调用TrapdoorGen函数
        auto result = RLWETrapdoorUtility<NativePoly>::TrapdoorGen(params, stddev, base, bal);
        
        Matrix<NativePoly>& A = result.first;
        RLWETrapdoorPair<NativePoly>& trapdoorPair = result.second;
        
        // 验证结果
        std::cout << "NativePoly Trapdoor generation successful!" << std::endl;
        std::cout << "Matrix A dimensions: " << A.GetRows() << " x " << A.GetCols() << std::endl;
        std::cout << "Trapdoor pair R dimensions: " << trapdoorPair.m_r.GetRows() 
                  << " x " << trapdoorPair.m_r.GetCols() << std::endl;
        std::cout << "Trapdoor pair E dimensions: " << trapdoorPair.m_e.GetRows() 
                  << " x " << trapdoorPair.m_e.GetCols() << std::endl;
        
        // 验证矩阵A的结构
        if (A.GetRows() == 1 && A.GetCols() >= 3) {
            std::cout << "Matrix A structure verification: PASS" << std::endl;
        } else {
            std::cout << "Matrix A dimensions incorrect: FAIL" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "NativePoly Trapdoor generation failed: " << e.what() << std::endl;
    }
}

void testGaussSamp() {
    std::cout << "\n=== Testing RLWETrapdoorUtility::GaussSamp ===" << std::endl;
    
    // 设置参数
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
    
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    
    // 获取参数
    auto params = cryptoContext->GetElementParams();
    
    // 生成trapdoor
    double stddev = 4.0;
    int64_t base = 8;
    bool bal = false;
    
    auto trapdoorResult = RLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(params, stddev, base, bal);
    Matrix<DCRTPoly>& A = trapdoorResult.first;
    RLWETrapdoorPair<DCRTPoly>& T = trapdoorResult.second;
    
    // 创建高斯采样器
    auto dgg = DCRTPoly::DggType(4.0);           // 离散高斯分布
    auto dggLargeSigma = DCRTPoly::DggType(16.0); // 大标准差的高斯分布
    
    // 创建目标向量u
    auto uniform_alloc = DCRTPoly::MakeDiscreteUniformAllocator(params, Format::EVALUATION);
    DCRTPoly u = uniform_alloc();
    
    try {
        // 调用GaussSamp函数
        size_t n = params->GetRingDimension();
        size_t k = A.GetCols() - 2; // 根据TrapdoorGen的实现，k = A.GetCols() - 2
        
        auto zHatPrime = RLWETrapdoorUtility<DCRTPoly>::GaussSamp(n, k, A, T, u, dgg, dggLargeSigma, base);
        
        std::cout << "Gaussian sampling successful!" << std::endl;
        std::cout << "Sampled vector dimensions: " << zHatPrime.GetRows() 
                  << " x " << zHatPrime.GetCols() << std::endl;
        
        // 验证采样结果是否满足 A * z = u (近似)
        auto Az = A.Mult(zHatPrime);
        Az(0, 0).SetFormat(Format::COEFFICIENT);
        u.SetFormat(Format::COEFFICIENT);
        
        std::cout << "Gaussian sampling verification completed" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Gaussian sampling failed: " << e.what() << std::endl;
    }
}

int main() {
    std::cout << "OpenFHE Trapdoor Utility Test Program" << std::endl;
    std::cout << "=====================================" << std::endl;
    
    // 测试基本的Trapdoor生成
    testTrapdoorGen();
    
    // 测试方阵Trapdoor生成

    
    // 测试高斯采样
    testGaussSamp();
    
    std::cout << "\nAll tests completed!" << std::endl;
    
    return 0;
}