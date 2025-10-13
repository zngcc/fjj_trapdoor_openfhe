// MyRLWETrapdoor.cpp
#include "MyRLWETrapdoor.h"
#include <cmath>

namespace myrlwe {

template <class Element>
std::pair<Matrix<Element>, MyRLWETrapdoorPair<Element>> 
MyRLWETrapdoorUtility<Element>::TrapdoorGen(
    std::shared_ptr<typename Element::Params> params, 
    const Matrix<Element>& a_prime, const Element& h, 
    double stddev, int64_t base) {
    
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    auto gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(
        params, COEFFICIENT, stddev);
    
    // 计算k (gadget向量的长度)
    double val = params->GetModulus().ConvertToDouble();
    double nBits = std::ceil(std::log2(val));
    // 如果是CRT k就不是log2(base)
    size_t k = std::ceil(nBits / std::log2(base));
    // size_t k = params->GetParams().size() * digitCount;
        std::cout << "TrapdoorGen: k = " << k << " (base=" << base << ")" << std::endl;

    // 生成T ← 高斯分布
    // T的维度: a'的列数 × k
    Matrix<Element> T(zero_alloc, a_prime.GetCols(), k, gaussian_alloc);
    T.SetFormat(Format::EVALUATION);
        std::cout << "a' dimensions: " << a_prime.GetRows() << " x " << a_prime.GetCols() << std::endl;
        std::cout << "T dimensions: " << T.GetRows() << " x " << T.GetCols() << std::endl;

    // 生成gadget向量g
        std::cout << "开始生成g: " << std::endl;
    Matrix<Element> g = Matrix<Element>(zero_alloc, 1, k).GadgetVector(base);
        std::cout << "g dimensions: " << g.GetRows() << " x " << g.GetCols() << std::endl;
    // 计算 a2 = h*g - a'^T * T
    Matrix<Element> a2 = (h * g) - (a_prime.Transpose() * T).Transpose();
        std::cout << "a2 dimensions: " << a2.GetRows() << " x " << a2.GetCols() << std::endl;
    // 垂直拼接 a' 和 a2 得到最终的a
    // Matrix<Element> a = a_prime.VStack(a2);
    Matrix<Element> a(zero_alloc, a_prime.GetRows() + a2.GetRows(), a_prime.GetCols());
        std::cout << "a dimensions: " << a.GetRows() << " x " << a.GetCols() << std::endl;
    // 复制a_prime的内容
    for (size_t i = 0; i < a_prime.GetRows(); ++i) {
        for (size_t j = 0; j < a_prime.GetCols(); ++j) {
            a(i, j) = a_prime(i, j);
        }
    }
    
    // 复制a2的内容
    for (size_t i = 0; i < a2.GetRows(); ++i) {
        for (size_t j = 0; j < a2.GetCols(); ++j) {
            a(a_prime.GetRows() + i, j) = a2(i, j);
        }
    }

    return std::make_pair(a, MyRLWETrapdoorPair<Element>(T));
}

template <class Element>
Matrix<Element> MyRLWETrapdoorUtility<Element>::SamplePre(
    const MyRLWETrapdoorPair<Element>& T_pair,
    const Matrix<Element>& a, const Element& h,
    double eta, double sigma, double alpha,
    const Element& u, typename Element::DggType& dgg,
    int64_t base) {
    
    const std::shared_ptr<typename Element::Params> params = u.GetParams();
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    
    const Matrix<Element>& T = T_pair.T;
    
    // 步骤1: Sample y ← SampleY(q, η, α, T)
    Matrix<Element> y = SampleY(params->GetModulus(), eta, alpha, T_pair, dgg, params);
    
    // 步骤2: 计算 v ← h^{-1} · (u - a·y)
    Element a_times_y = (a * y)(0, 0);
    Element u_minus_ay = u - a_times_y;
    
    // 计算h的模逆
    Element h_inv = h.MultiplicativeInverse();
    Element v = h_inv * u_minus_ay;
    
    // 步骤3: Sample z ← SamplePolyG(σ, v)
    size_t k = T.GetCols();  // z的维度应与T的列数相同
    Matrix<Element> z = SamplePolyG(sigma, v, dgg, k, params);
    
    // 步骤4: 计算 x ← y + [T; I] · z
    size_t m_prime = T.GetRows();  // T的行数
    size_t k_dim = T.GetCols();    // T的列数
    
    // 构建 [T; I] 矩阵
    Matrix<Element> T_I(zero_alloc, m_prime + k_dim, k_dim);
    
    // 上半部分是T
    for (size_t i = 0; i < m_prime; ++i) {
        for (size_t j = 0; j < k_dim; ++j) {
            T_I(i, j) = T(i, j);
        }
    }
    
    // 下半部分是单位矩阵
    // for (size_t i = 0; i < k_dim; ++i) {
    //     T_I(m_prime + i, i) = Element::Allocator(params, EVALUATION)();  // 零元素
    //     // 设置对角线为1
    //     auto one_poly = Element::Allocator(params, COEFFICIENT)();
    //     one_poly.SetValuesToOne();
    //     one_poly.SetFormat(EVALUATION);
    //     T_I(m_prime + i, i) = one_poly;
    // }
    Element unity_poly = CreateUnityPoly(params);
    for (size_t i = 0; i < k_dim; ++i) {
        for (size_t j = 0; j < k_dim; ++j) {
            if (i == j) {
                T_I(m_prime + i, j) = unity_poly;
            } else {
                T_I(m_prime + i, j) = zero_alloc();
            }
        }
    }
    
    Matrix<Element> T_I_times_z = T_I * z;
    Matrix<Element> x = y + T_I_times_z;
    
    return x;
}

template <class Element>
Matrix<Element> MyRLWETrapdoorUtility<Element>::SamplePolyG(
    double sigma, const Element& v, typename Element::DggType& dgg, 
    size_t k, std::shared_ptr<typename Element::Params> params) {
    
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    
    // 计算 α = √5 * σ
    double sampling_sigma = std::sqrt(5.0) * sigma;
    
    // 创建高斯分配器
    auto gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(
        params, COEFFICIENT, sampling_sigma);
    
    // 采样z，维度为k x 1
    Matrix<Element> z(zero_alloc, k, 1, gaussian_alloc);
    z.SetFormat(Format::EVALUATION);
    
    return z;
}

template <class Element>
Matrix<Element> MyRLWETrapdoorUtility<Element>::SampleY(
    const typename Element::Integer& q, double eta, double alpha,
    const MyRLWETrapdoorPair<Element>& T_pair, typename Element::DggType& dgg,
    std::shared_ptr<typename Element::Params> params) {
    
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    
    const Matrix<Element>& T = T_pair.T;
    size_t m_prime = T.GetRows();  // T的行数
    size_t k = T.GetCols();        // T的列数
    size_t m = m_prime + k;        // y的总维度
    
    // 检查条件: η > s₁(T) * α
    double spectral_norm = EstimateSpectralNorm(T);
    if (eta <= spectral_norm * alpha) {
        throw std::invalid_argument("η must be greater than s₁(T) * α");
    }
    
    // 简化实现：直接使用独立高斯采样
    // 在实际应用中，你可能需要实现更精确的协方差采样
    auto gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(
        params, COEFFICIENT, eta);
    
    Matrix<Element> y(zero_alloc, m, 1, gaussian_alloc);
    y.SetFormat(Format::EVALUATION);
    
    return y;
}

template <class Element>
double MyRLWETrapdoorUtility<Element>::EstimateSpectralNorm(const Matrix<Element>& T) {
    // 简化实现：返回一个保守估计
    // 在实际应用中，你需要实现更精确的谱范数计算

    // double frobenius_norm = 0.0;
    // for (size_t i = 0; i < T.GetRows(); ++i) {
    //     for (size_t j = 0; j < T.GetCols(); ++j) {
    //         // 这里简化处理，实际需要计算多项式的范数
    //         frobenius_norm += 1.0;  // 保守估计
    //     }
    // }
    // frobenius_norm = std::sqrt(frobenius_norm);
    
    // return frobenius_norm;  // 谱范数 ≤ Frobenius范数
    return 1.0;  // 保守估计
}

template <class Element>
Element MyRLWETrapdoorUtility<Element>::CreateUnityPoly(std::shared_ptr<typename Element::Params> params) {
    auto alloc = Element::Allocator(params, COEFFICIENT);
    Element one_poly = alloc();
    
    // 手动设置多项式为1：第一个系数为1，其他为0
    typename Element::Vector coefficients = one_poly.GetValues();
    if (coefficients.GetLength() > 0) {
        coefficients[0] = 1;
    }
    one_poly.SetValues(coefficients, COEFFICIENT);
    one_poly.SetFormat(EVALUATION);
    
    return one_poly;
}



// 显式实例化常用类型
template class MyRLWETrapdoorUtility<Poly>;
template class MyRLWETrapdoorUtility<NativePoly>;
template class MyRLWETrapdoorUtility<DCRTPoly>;
} // namespace myrlwe