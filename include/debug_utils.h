// debug_utils.h
#ifndef DEBUG_UTILS_H
#define DEBUG_UTILS_H

#include <iostream>
#include <openfhe.h>

namespace debug {

using namespace lbcrypto;

// ==================== 函数1：识别并打印类型信息 ====================

/**
 * 识别输入类型并打印相关信息
 * 如果是矩阵：打印矩阵大小和元素类型
 * 如果是单个元素：打印元素类型
 */
template<typename T>
inline void printTypeInfo(const T& obj, const std::string& label = "") {
    if constexpr (std::is_same_v<T, Matrix<DCRTPoly>>) {
        // 处理DCRTPoly矩阵
        std::cout <<  label << ":" << obj.GetRows() << " x " << obj.GetCols();
        
        // 打印矩阵中一个元素的详细信息
        if (obj.GetRows() > 0 && obj.GetCols() > 0) {
            const DCRTPoly& elem = obj(0, 0);
            // std::cout << " type: DCRTPoly" ;
            std::cout << " ring dim: " << elem.GetRingDimension();
            std::cout << " CRT moduli count: " << elem.GetNumOfElements();
            std::cout << " format: " << (elem.GetFormat() == Format::COEFFICIENT ? "COEFFICIENT" : "EVALUATION") << std::endl;
        }
    }
    else if constexpr (std::is_same_v<T, DCRTPoly>) {
        // 处理单个DCRTPoly
        std::cout <<  label << ":" << "1 x 1";
        // std::cout << " type: DCRTPoly" ;
        std::cout << " ring dim: " << obj.GetRingDimension();
        std::cout << " CRT moduli count: " << obj.GetNumOfElements();
        std::cout << " format: " << (obj.GetFormat() == Format::COEFFICIENT ? "COEFFICIENT" : "EVALUATION") << std::endl;
    }
    else {
        // 其他类型
        std::cout << label << " 1 x 1; type: " << typeid(T).name() << std::endl;
    }
    std::cout << std::endl;
}

// ==================== 函数2：转换并打印DCRTPoly为普通多项式 ====================

/**
 * 将DCRTPoly转换为普通的NativePoly（使用第一个模数）
 * 用于调试查看多项式系数
 */
inline NativePoly convertToNativePoly(const DCRTPoly& dcrtPoly, size_t model_id = 2) {
    // 复制一份以确保不改变原多项式格式
    DCRTPoly temp = dcrtPoly;
    
    // 确保在系数格式下
    if (temp.GetFormat() != Format::COEFFICIENT) {
        temp.SetFormat(Format::COEFFICIENT);
    }
    
    // 检查是否有CRT分量
    if (temp.GetNumOfElements() == 0) {
        throw std::runtime_error("DCRTPoly没有CRT分量");
    }
    
    // 获取第一个模数的参数
    auto firstParams = temp.GetParams()->GetParams()[model_id];
    
    // 获取第一个CRT分量的系数
    auto firstElement = temp.GetElementAtIndex(model_id);
    
    // 创建NativePoly
    NativePoly nativePoly(firstParams, Format::COEFFICIENT, false);
    
    // 获取系数并设置到NativePoly
    auto coefficients = firstElement.GetValues();
    nativePoly.SetValues(coefficients, Format::COEFFICIENT);
    
    return nativePoly;
}

/**
 * 转换并打印单个DCRTPoly
 * @param poly 要转换的DCRTPoly
 * @param label 标签（可选）
 * @param numCoeffs 要打印的系数数量（默认为10个）
 */
inline void printAsNormalPoly(const DCRTPoly& poly, 
                      const std::string& label = "DCRTPoly", 
                      size_t numCoeffs = 10) {
    std::cout << label ;
    
    try {
        // 获取模数信息
        auto params = poly.GetParams();
        if (!params) {
            std::cout << "错误: 无法获取参数" << std::endl;
            return;
        }
        
        // 输出模数信息
        size_t numModuli = params->GetParams().size();
        std::cout << " "<< numModuli << "个模数: ";
        
        // 显示所有模数
        for (size_t i = 0; i < numModuli; ++i) {
            auto modParam = params->GetParams()[i];
            auto modulus = modParam->GetModulus();
            std::cout <<  modulus;
            if (i < numModuli - 1) std::cout << ", ";
        }
        
        // 计算总模数（所有模数的乘积）
        if (numModuli > 0) {
            auto totalModulus = params->GetParams()[0]->GetModulus();
            for (size_t i = 1; i < numModuli; ++i) {
                totalModulus = totalModulus * params->GetParams()[i]->GetModulus();
            }
            std::cout << " 总模数: " << totalModulus;
            NativePoly nativePoly =  convertToNativePoly(poly, 0);
            std::cout << " 多项式环维度: " << nativePoly.GetRingDimension();
            std::cout << " 格式: " << (nativePoly.GetFormat() == Format::COEFFICIENT ? "COEFFICIENT" : "EVALUATION") << std::endl;
                        // 打印系数
            for (size_t i = 0; i < numModuli; ++i) {
                nativePoly = convertToNativePoly(poly, i);
                std::cout << "Module " << i << " :";
                auto coeffs = nativePoly.GetValues();
                size_t limit = std::min(numCoeffs, coeffs.GetLength());
                
                for (size_t i = 0; i < limit; ++i) {
                    std::cout << " [" << i << "]:" << coeffs[i];
                }
                
                if (coeffs.GetLength() >= numCoeffs) {
                    std::cout << "  ... 余 " << (coeffs.GetLength() - numCoeffs) << " 个" << std::endl;
                }
            }
        }

    }
    catch (const std::exception& e) {
        std::cout << "转换失败: " << e.what() << std::endl;
    }
    
    // std::cout << std::endl;
}
/**
 * 转换并打印DCRTPoly矩阵
 * @param matrix 要转换的Matrix<DCRTPoly>
 * @param label 标签（可选）
 * @param numElements 要转换和打印的元素数量（默认为前3个）
 * @param numCoeffs 每个元素要打印的系数数量（默认为5个）
 */
template<typename T>
inline void printAsNormalPoly(const T& matrix,
                      const std::string& label = "DCRTPoly矩阵", 
                      size_t numElements = 3,
                      size_t numCoeffs = 5) {
    if constexpr (std::is_same_v<T, Matrix<DCRTPoly>>) {
        std::cout << label << " " << matrix.GetRows() << " x " << matrix.GetCols() << ":" << std::endl;
        
        // 计算要打印的元素数量
        size_t totalElements = matrix.GetRows() * matrix.GetCols();
        size_t elementsToPrint = std::min(numElements, totalElements);
        
        // std::cout << "显示前" << elementsToPrint << "个元素的转换结果:" << std::endl;
        
        // 遍历并转换前几个元素
        size_t count = 0;
        for (size_t i = 0; i < matrix.GetRows() && count < elementsToPrint; ++i) {
            for (size_t j = 0; j < matrix.GetCols() && count < elementsToPrint; ++j) {
                std::cout << "--- 元素 [" << i << "," << j << "] ---" << std::endl;
                
                try {
                    // 转换为NativePoly并打印
                    printAsNormalPoly(matrix(i, j), "", numCoeffs);
                }
                catch (const std::exception& e) {
                    std::cout << "转换失败: " << e.what() << std::endl;
                }
                
                count++;
            }
        }
    }
    else if constexpr (std::is_same_v<T, DCRTPoly>) {
        std::cout << label << " (单元素):" << std::endl;
        
        try {
            printAsNormalPoly(matrix, "", numCoeffs);
        }
        catch (const std::exception& e) {
            std::cout << "转换失败: " << e.what() << std::endl;
        }
    }
    else {
        // 其他类型
        std::cout << label << " 1 x 1; type: " << typeid(T).name() << std::endl;
    }
    std::cout << std::endl;
}

} // namespace debug

#endif // DEBUG_UTILS_H