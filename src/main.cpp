#include <iostream>
#include "util.h"
#include <openfhe.h>
// main.cpp
#include "MyRLWETrapdoor.h"
#include "debug_utils.h"
using namespace lbcrypto;
using namespace myrlwe;


#include <vector>
#include <complex>


template <class Element>
void demonstrateMyRLWEDCRTPoly() {
    // 1. 初始化OpenFHE上下文
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 15;
    uint32_t batchSize = 4;
    // SecurityLevel securityLevel = HEStd_128_classic;
    SecurityLevel securityLevel = HEStd_NotSet; //其他安全等级对toy参数都会抛出error
    // 创建加密参数
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleFactorBits);
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(securityLevel);
    parameters.SetRingDim(256);

    
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
    auto gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(params, Format::COEFFICIENT, 4.0);
    auto uniform_alloc = Element::MakeDiscreteUniformAllocator(params, EVALUATION);
        std::cout << "开始在外面生成k"<< std::endl;
    NativeInteger q = params->GetParams()[0]->GetModulus();
    size_t m = params->GetRingDimension();
    size_t digitCount = static_cast<size_t>(ceil(log2(q.ConvertToDouble()) / log2(2)));

    size_t k = params->GetParams().size() * digitCount;
        std::cout << "k = " << k << " m = " << m << std::endl;

    Matrix<Element> a_prime(zero_alloc, m-k, 1);
        std::cout << "创建a'空间" << std::endl;
    for (size_t j = 0; j < m-k; ++j) {
        a_prime(j, 0) = uniform_alloc();
    }
        std::cout << "生成了a'" << std::endl;

    // 5. 创建h元素
    DCRTPoly h(params, Format::EVALUATION);
    h = 2;  // 设置为常数
        debug::printAsNormalPoly(h, "h", 10);
    Element h_inv = h.MultiplicativeInverse();
        debug::printAsNormalPoly(h_inv, "h_inv", 10);


    Element hhinv = h * h_inv;
        debug::printAsNormalPoly(hhinv, "hhinv", 10);
    // 6. 使用我们的Trapdoor生成算法
    double sigma = 3.2;
    auto trapdoor_result = MyRLWETrapdoorUtility<Element>::TrapdoorGen(
        params, a_prime, h, sigma);
    
    Matrix<Element> A = trapdoor_result.first;
    Matrix<Element> T = trapdoor_result.second;
    // a = [a'^T| hg-a'^TT]^T
        debug::printTypeInfo(A, "A");
        debug::printTypeInfo(T, "T");
    // x ← y + [T; I] · (h^-1 * (u - a'^Ty))
    std::cout << "Trapdoor generation completed!" << std::endl;
    
    // 7. 使用SamplePre算法
    Element u = uniform_alloc();  // 目标syndrome
    double eta = 8.0,  alpha = 1.0;
    
    std::cout << "Begin Sampling Preimage" << std::endl;
    Matrix<Element> x = MyRLWETrapdoorUtility<Element>::SamplePre(
        T, A, h, eta, sigma, alpha, u);
    
    std::cout << "SamplePre completed!" << std::endl;
    
    // 8. 验证: 检查 A·x ≈ h·u
    debug::printTypeInfo(A, "A");
    debug::printTypeInfo(x, "x");
    Matrix<Element> Ax = A.Transpose() * x;
        debug::printTypeInfo(Ax, "Ax");

        debug::printTypeInfo(h, "h");
        debug::printTypeInfo(u, "u");

    // 计算误差
    Element error = Ax(0, 0) - u;
        debug::printAsNormalPoly(error, "error");
    // error.SetFormat(COEFFICIENT);
    //     debug::printAsNormalPoly(error, "error");

    std::cout << "Verification completed!" << std::endl;
    // 可以进一步分析error的大小来验证正确性
}

template <class Element>
void MyRLWEDCRTPolyEncDec() {
    // 1. 初始化OpenFHE上下文
    uint32_t multDepth = 1;
    uint32_t scaleFactorBits = 15;

    uint32_t batchSize = 4;
    // SecurityLevel securityLevel = HEStd_128_classic;
    SecurityLevel securityLevel = HEStd_NotSet; //其他安全等级对toy参数都会抛出error
    // 创建加密参数
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleFactorBits);
    parameters.SetBatchSize(batchSize);
    parameters.SetSecurityLevel(securityLevel);
    parameters.SetRingDim(256);

    
    CryptoContext<Element> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    
    // 3. 获取环参数
    auto params = cc->GetElementParams();
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    auto gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(params, Format::COEFFICIENT, 3.2);
    auto uniform_alloc = Element::MakeDiscreteUniformAllocator(params, EVALUATION);
        std::cout << "开始在外面生成k"<< std::endl;
    NativeInteger q = params->GetParams()[0]->GetModulus();
    size_t m = params->GetRingDimension();
    size_t digitCount = static_cast<size_t>(ceil(log2(q.ConvertToDouble()) / log2(2)));

    size_t k = params->GetParams().size() * digitCount;
        std::cout << "k = " << k << " m = " << m << std::endl;
            // ==================== Setup 算法 ====================
    std::cout << "\n=== 执行 Setup 算法 ===" << std::endl;
    Element onepoly(params, Format::EVALUATION);
    onepoly = 1;  // 设置为常数1
    Matrix<Element> a_root(zero_alloc, m-k, 1);
        std::cout << "创建a'空间" << std::endl;
        a_root(0,0) = onepoly;
    for (size_t j = 1; j < m-k; ++j) {
        a_root(j, 0) = uniform_alloc();
    }
        std::cout << "生成了a'" << std::endl;
        debug::printAsNormalPoly(a_root, "a_root", 10);

    double sigma = 3.2;
    Matrix<DCRTPoly> s_root = MyRLWETrapdoorUtility<DCRTPoly>::GenSecretKey(m - k, sigma, params);
        debug::printAsNormalPoly(s_root, "root secret key");

    auto as = a_root.Transpose() * s_root;  // 目标syndrome
    Element u = as(0,0);

    std::cout << "Create s " << std::endl;
    Matrix<Element> s(zero_alloc, m, 1);
    for (size_t i = 0; i < m-k; ++i) {
        s(i, 0) = s_root(i, 0);
    }
        debug::printAsNormalPoly(s, "root secret key"); //75
    ////////////////////////////////////////////////////
    // evalute secret key
    // std::cout << "Create evalute secret key " << std::endl;
    // size_t evk_size = (m-k) * (m-k);
    // std::vector<std::pair<Element, Matrix<Element>>> evk_par(evk_size);
    // Matrix<Element> a_prime_evk(zero_alloc, m, 1);
    // for (size_t i = 0; i < m; ++i) {
    //     a_prime_evk(i, 0) = uniform_alloc();
    // }
    // for (size_t idx = 0; idx < evk_size; ++idx) {
    //     size_t i = idx / (m-k);
    //     size_t j = idx % (m-k);
        
    //     // 随机选择 P_{i,j}
    //     int64_t P_ij = (std::rand() % 100) + 1;  // 1-100之间的随机整数
        
    //     // 获取 s_i 和 s_j
    //     Element s_i = s_root(i, 0);
    //     Element s_j = s_root(j, 0);
        
    //     // 计算 aux_{i,j} = P_{i,j} * s_i * s_j
    //     Element aux_ij = s_i * s_j;
    //     // 这里简化处理，实际需要乘以P_{i,j}，但由于环运算，需要将P转换为多项式
        
    //     // 采样错误 e_{i,j}'
    //     Element e_ij_prime = gaussian_alloc();
    //     e_ij_prime.SetFormat(EVALUATION);
        
    //     // 计算 b'_{i,j} = -a'·s + e_{i,j}' + aux_{i,j}
    //     auto a_prime_dot_s = a_prime_evk.Transpose() * s;
    //     Element b_ij_prime = -a_prime_dot_s(0, 0) + e_ij_prime + aux_ij;
        
    //     // 存储评估密钥
    //     evk_par[idx] = std::make_pair(b_ij_prime, a_prime_evk);
    // }

    ////////////////////////////////////////////////////
    
    // 5. 创建h元素
    DCRTPoly h(params, Format::EVALUATION);
    h = 2;  // 设置为常数
        debug::printAsNormalPoly(h, "h", 10);
    Element h_inv = h.MultiplicativeInverse();
        debug::printAsNormalPoly(h_inv, "h_inv", 10);


    Element hhinv = h * h_inv;
        debug::printAsNormalPoly(hhinv, "hhinv", 10);
    // 6. 使用我们的Trapdoor生成算法

    auto trapdoor_result = MyRLWETrapdoorUtility<Element>::TrapdoorGen(
        params, a_root, h, sigma);
    
    Matrix<Element> A = trapdoor_result.first;
    Matrix<Element> T = trapdoor_result.second;
    // a = [a'^T| hg-a'^TT]^T
        debug::printTypeInfo(A, "A");
        debug::printTypeInfo(T, "T");
    // x ← y + [T; I] · (h^-1 * (u - a'^Ty))
    std::cout << "Trapdoor generation completed!" << std::endl;
    
    // 7. 使用SamplePre算法

    double eta = 8.0, alpha = 1.0;
    
    
    std::cout << "Begin Sampling Preimage" << std::endl;
    Matrix<Element> x = MyRLWETrapdoorUtility<Element>::SamplePre(
        T, A, h, eta, sigma, alpha, u);
    
    std::cout << "SamplePre completed!" << std::endl;
    
    // 8. 验证: 检查 A·x ≈ u
    debug::printTypeInfo(A, "A");
    debug::printTypeInfo(x, "x");
    Matrix<Element> Ax = A.Transpose() * x;
        debug::printTypeInfo(Ax, "Ax");

        debug::printTypeInfo(h, "h");
        debug::printTypeInfo(u, "u");

    // 计算误差
    Element error = Ax(0, 0) - u;
        debug::printAsNormalPoly(error, "error");
    // error.SetFormat(COEFFICIENT);
    //     debug::printAsNormalPoly(error, "error");

    std::cout << "Verification completed!" << std::endl;
    

    // 5. 生成子用户的评估密钥 evk_child (与evk_par类似，但使用新的a')
    // std::vector<std::pair<Element, Matrix<Element>>> evk_child(evk_size);
    // Matrix<Element> a_prime_child(zero_alloc, m, 1);
    // for (size_t i = 0; i < m; ++i) {
    //     a_prime_child(i, 0) = uniform_alloc();
    // }
    //     for (size_t idx = 0; idx < evk_size; ++idx) {
    //     size_t i = idx / (m-k);
    //     size_t j = idx % (m-k);
        
    //     int64_t P_ij = (std::rand() % 100) + 1;
    //     Element s_i = s_root(i, 0);
    //     Element s_j = s_root(j, 0);
    //     Element aux_ij = s_i * s_j;
        
    //     Element e_ij_prime = gaussian_alloc();
    //     e_ij_prime.SetFormat(EVALUATION);
        
    //     auto a_prime_dot_s = a_prime_child.Transpose() * s;
    //     Element b_ij_prime = -a_prime_dot_s(0, 0) + e_ij_prime + aux_ij;
        
    //     evk_child[idx] = std::make_pair(b_ij_prime, a_prime_child);
    // }

    // ==================== Encrypt 算法 ====================
    std::cout << "\n=== 执行 Encrypt 算法 ===" << std::endl;
    auto encodingParams = cc->GetEncodingParams();
    // BatchEncoding batchEncoder(cc);
    
    // 创建明文向量
    // std::vector<double> plaintextVector = {1.0, 2.0, 3.0, 4.0};

    // DCRTPoly h(params, Format::EVALUATION);
    DCRTPoly mu(params, Format::EVALUATION);
    // mu = 5;
    std::vector<int64_t> coeff = {5};
    mu.SetCoefficients(coeff);
    debug::printAsNormalPoly(mu, "mu");
    
        // Element scaleFactor(params, Format::EVALUATION);
    float scaleFactor = pow(2.0, scaleFactorBits);  // CKKS缩放因子
    auto scaled_mu = mu * scaleFactor;
    debug::printAsNormalPoly(scaled_mu, "scaled mu");
    // 采样随机多项式 p
    DCRTPoly p(params, Format::EVALUATION);
    // p = 2;  // 设置为常数
    p = uniform_alloc();
    
    // 采样错误向量
        double small_sigma = 0.5;  // 更小的标准差
    auto small_gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(
        params, Format::COEFFICIENT, small_sigma);
    Matrix<Element> e0(zero_alloc, m-k, 1, small_gaussian_alloc);
    e0.SetFormat(EVALUATION);
    Matrix<Element> e1(zero_alloc, k, 1, small_gaussian_alloc);
    e1.SetFormat(EVALUATION);
    Matrix<Element> e(zero_alloc, m, 1);
    for (size_t i = 0; i < m-k; ++i) {
        e(i, 0) = e0(i, 0);
    }
    for (size_t i = 0; i < k; ++i) {
        e(m-k + i, 0) = e1(i, 0);
    }
    // 采样错误 e2
    Element e2 = small_gaussian_alloc();
    e2.SetFormat(EVALUATION);
        // std::cout << "get e" << std::endl;
        debug::printAsNormalPoly(e2, "e2");
                debug::printAsNormalPoly((e2)*scaleFactor.MultiplicativeInverse(), "e2 /scale");
                debug::printAsNormalPoly((e2 * scaleFactor)*scaleFactor.MultiplicativeInverse(), "(scaled e2) /scale");
        debug::printAsNormalPoly((e2 + scaled_mu)*scaleFactor.MultiplicativeInverse(), "e2 + scaled mu /scale");
        debug::printAsNormalPoly(mu, "mu");
                debug::printAsNormalPoly(scaled_mu*scaleFactor.MultiplicativeInverse(), "scaled mu /scale");
    // 计算密文
    // c0 = u·p + e2 + μ
    Element c0 = u * p + e2 + scaled_mu;
        // std::cout << "get c0" << std::endl;
        debug::printAsNormalPoly(c0, "c0");
    // c1 = -A·p + e
        debug::printAsNormalPoly(A, "A");
        debug::printAsNormalPoly(p, "p");
        debug::printAsNormalPoly(e, "e");
    Matrix<Element> A_p = A * p;

    Matrix<Element> c1 = e - A_p;
            debug::printAsNormalPoly(c1, "c1");
    std::cout << "加密完成，生成密文 (c0, c1)" << std::endl;

        // ==================== ChildDec 算法 (CKKS版本) ====================
    std::cout << "\n=== 执行 CKKS ChildDec 算法 ===" << std::endl;
    
    // 计算 m_raw = c0 + c1·x
    Matrix<Element> c1_dot_x = c1.Transpose() * x;
        // debug::printAsNormalPoly(c1, "c1");
                // debug::printAsNormalPoly(x, "x");
    Element m_raw = c0 + c1_dot_x(0, 0);
        debug::printAsNormalPoly(m_raw, "m_raw");
        debug::printAsNormalPoly(m_raw-scaled_mu, "m_raw - scaled mu");

}
int main() {
    std::cout << "trapdooribe fjj" << std::endl;
    std::cout << "================" << std::endl;
    // std::cout << "Testing with Poly..." << std::endl;
    std::cout << "Testing with DCRTPoly..." << std::endl;
    // demonstrateMyRLWEDCRTPoly<DCRTPoly>();
    MyRLWEDCRTPolyEncDec<DCRTPoly>();
    return 0;
}
