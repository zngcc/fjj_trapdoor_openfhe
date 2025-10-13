#include "util.h"
#include <openfhe.h>

using namespace lbcrypto;

void CKKSBasicOperations() {
    std::cout << "=== OpenFHE CKKS 基础操作演示 ===" << std::endl;
    
    // 设置参数
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(4);
    parameters.SetScalingModSize(50);
    parameters.SetBatchSize(8);  // SIMD 槽位数量
    
    // 创建加密上下文
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    std::cout << "1. 加密上下文创建成功" << std::endl;
    
    // 启用功能
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    std::cout << "2. 所需功能已启用" << std::endl;
    
    // 生成密钥对
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
    std::cout << "3. 密钥对生成成功" << std::endl;
    
    // 准备明文数据 - 使用SIMD编码
    std::vector<double> plaintext1 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<double> plaintext2 = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    
    std::cout << "\n原始明文向量1: ";
    print_vector(plaintext1);
    std::cout << "原始明文向量2: ";
    print_vector(plaintext2);
    
    // 编码明文
    Plaintext pt1 = cryptoContext->MakeCKKSPackedPlaintext(plaintext1);
    Plaintext pt2 = cryptoContext->MakeCKKSPackedPlaintext(plaintext2);
    std::cout << "4. 明文编码成功 (SIMD)" << std::endl;
    
    // 加密
    Ciphertext<DCRTPoly> ct1 = cryptoContext->Encrypt(keyPair.publicKey, pt1);
    Ciphertext<DCRTPoly> ct2 = cryptoContext->Encrypt(keyPair.publicKey, pt2);
    std::cout << "5. 加密成功" << std::endl;
    
    Timer timer;
    
    // 同态加法操作
    std::cout << "\n=== 同态加法操作 ===" << std::endl;
    timer.start();
    auto ctAdd = cryptoContext->EvalAdd(ct1, ct2);
    double addTime = timer.stop();
    
    // 解密加法结果
    Plaintext ptAdd;
    cryptoContext->Decrypt(keyPair.secretKey, ctAdd, &ptAdd);
    ptAdd->SetLength(plaintext1.size());
    std::vector<double> resultAdd = ptAdd->GetRealPackedValue();
    
    std::cout << "加法结果: ";
    print_vector(resultAdd);
    std::cout << "加法耗时: " << addTime * 1000 << " ms" << std::endl;
    
    // 同态乘法操作
    std::cout << "\n=== 同态乘法操作 ===" << std::endl;
    timer.start();
    auto ctMult = cryptoContext->EvalMult(ct1, ct2);
    double multTime = timer.stop();
    
    // 解密乘法结果
    Plaintext ptMult;
    cryptoContext->Decrypt(keyPair.secretKey, ctMult, &ptMult);
    ptMult->SetLength(plaintext1.size());
    std::vector<double> resultMult = ptMult->GetRealPackedValue();
    
    std::cout << "乘法结果: ";
    print_vector(resultMult);
    std::cout << "乘法耗时: " << multTime * 1000 << " ms" << std::endl;
    
    // 验证结果正确性
    std::cout << "\n=== 结果验证 ===" << std::endl;
    
    // 计算预期结果
    std::vector<double> expectedAdd(plaintext1.size());
    std::vector<double> expectedMult(plaintext1.size());
    for (size_t i = 0; i < plaintext1.size(); ++i) {
        expectedAdd[i] = plaintext1[i] + plaintext2[i];
        expectedMult[i] = plaintext1[i] * plaintext2[i];
    }
    
    std::cout << "预期加法结果: ";
    print_vector(expectedAdd);
    std::cout << "预期乘法结果: ";
    print_vector(expectedMult);
    
    // 计算误差
    double addError = calculate_max_error(expectedAdd, resultAdd);
    double multError = calculate_max_error(expectedMult, resultMult);
    
    std::cout << "加法平均误差: " << addError << std::endl;
    std::cout << "乘法平均误差: " << multError << std::endl;
}

void CKKSWithBootstrap() {
    std::cout << "\n\n=== OpenFHE CKKS 包含Bootstrap的操作演示 ===" << std::endl;
    
    // 设置支持Bootstrap的参数 - 需要更完整的配置
    CCParams<CryptoContextCKKSRNS> parameters;
    
    // 1. 设置基础参数
    parameters.SetMultiplicativeDepth(5);
    parameters.SetScalingModSize(59);
    parameters.SetBatchSize(8);
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    
    // 2. 专门设置Bootstrapping参数
    parameters.SetNumLargeDigits(3);
    
    // 3. 确保设置Bootstrapping所需的环维度（根据您的硬件调整）
    // 常见的环维度：4096, 8192, 16384
    parameters.SetRingDim(1 << 14); // 8192
    
    // 4. 创建加密上下文
    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    
    // 5. 启用功能 - 注意顺序
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);
    cryptoContext->Enable(FHE); // 最后启用FHE
    
    std::cout << "1. 支持Bootstrap的加密上下文创建成功" << std::endl;
    
    // 6. 生成密钥对
    KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);
        // 8. 现在设置Bootstrapping
    cryptoContext->EvalBootstrapSetup();
    cryptoContext->EvalBootstrapKeyGen(keyPair.secretKey, parameters.GetBatchSize());
    // 7. 在调用EvalBootstrapSetup之前生成旋转密钥
    std::vector<int32_t> indices = {1, 2, 3, 4, 5, 6, 7}; // 根据需要调整
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, indices);
    

    
    std::cout << "2. Bootstrap密钥生成成功" << std::endl;
    
    // 测试数据
    std::vector<double> plaintext = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    Plaintext pt = cryptoContext->MakeCKKSPackedPlaintext(plaintext);
    auto ct = cryptoContext->Encrypt(keyPair.publicKey, pt);
    
    std::cout << "原始数据: ";
    print_vector(plaintext);
    
    // 进行多次乘法（消耗层级）
    Timer timer;
    auto ctResult = ct;
    for (int i = 0; i < 3; ++i) {
        ctResult = cryptoContext->EvalMult(ctResult, ct);
    }
    
    std::cout << "3. 完成多次乘法操作" << std::endl;
    
    // 如果需要，这里可以执行Bootstrap来刷新层级
    // auto ctBootstrapped = cryptoContext->EvalBootstrap(ctResult);
    
    // 解密结果
    Plaintext ptResult;
    cryptoContext->Decrypt(keyPair.secretKey, ctResult, &ptResult);
    ptResult->SetLength(plaintext.size());
    auto result = ptResult->GetRealPackedValue();
    
    std::cout << "多次乘法后结果: ";
    print_vector(result);
    
    // 计算预期结果
    std::vector<double> expected(plaintext.size());
    for (size_t i = 0; i < plaintext.size(); ++i) {
        expected[i] = plaintext[i] * plaintext[i] * plaintext[i] * plaintext[i];
    }
    
    std::cout << "预期结果: ";
    print_vector(expected);
    std::cout << "实际误差: " << calculate_max_error(expected, result) << std::endl;
}

int main() {
    try {
        // 运行基础操作演示
        CKKSBasicOperations();
        
        // 运行包含Bootstrap的演示
        CKKSWithBootstrap();
        
        std::cout << "\n=== 所有测试完成 ===" << std::endl;
        
    } catch (std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

