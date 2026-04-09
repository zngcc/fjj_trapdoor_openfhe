// MyRLWETrapdoor.cpp
#include "MyRLWETrapdoor.h"
#include <cmath>
#include "debug_utils.h"
namespace myrlwe {

template <>
std::pair<Matrix<DCRTPoly>, Matrix<DCRTPoly>> MyRLWETrapdoorUtility<DCRTPoly>::TrapdoorGen(
    std::shared_ptr<typename DCRTPoly::Params> params, 
    const Matrix<DCRTPoly>& a_prime, const DCRTPoly& h, 
    double sigma, int64_t base) {
    
    auto zero_alloc = DCRTPoly::Allocator(params, Format::EVALUATION);
    auto gaussian_alloc = DCRTPoly::MakeDiscreteGaussianCoefficientAllocator(
        params, Format::COEFFICIENT, sigma);
    
    // 计算k (gadget向量的长度)
    NativeInteger q = params->GetParams()[0]->GetModulus();
        std::cout << "q = " << q << " " << params << std::endl;
    size_t digitCount = static_cast<size_t>(ceil(log2(q.ConvertToDouble()) / log2(base)));

    size_t k = params->GetParams().size() * digitCount;

    // 如果是CRT k就不是log2(base)
        std::cout << "TrapdoorGen: k = " << k << " (base=" << base << ")" << std::endl;

    // 生成T ← 高斯分布
    // T的维度:  (m-k)× k
        debug::printAsNormalPoly(a_prime, "a'");
    Matrix<DCRTPoly> T(zero_alloc, a_prime.GetRows(), k, gaussian_alloc);
    T.SetFormat(Format::EVALUATION);
        debug::printAsNormalPoly(T, "T");

        Matrix<DCRTPoly> T_row0(zero_alloc, 1, T.GetCols());
        T_row0.SetFormat(Format::EVALUATION);
        for (size_t j = 0; j < T.GetCols(); ++j) {
            T_row0(0, j) = T(0, j);
        }
        debug::printAsNormalPoly(T_row0, "T[0]");


        Matrix<DCRTPoly> T_row1(zero_alloc, 1, T.GetCols());
        T_row1.SetFormat(Format::EVALUATION);
        for (size_t j = 0; j < T.GetCols(); ++j) {
            T_row1(0, j) = T(1, j);
        }
        debug::printAsNormalPoly(T_row1, "T[1]");
///////////////////////////////////////////
    const auto& moduli = params->GetParams();
    size_t numModuli = moduli.size();

        // Matrix<DCRTPoly> g(zero_alloc, 1, k);
    
    // 为每个模数生成gadget元素
    // 生成gadget向量g (维度k)
        std::cout << "开始生成g: " << std::endl;
        Matrix<DCRTPoly> gadget = GenGadget(k, params);
        
    // DCRTPoly onepoly(params, Format::EVALUATION);
    // onepoly = 1;  // 设置为常数
    // DCRTPoly twopoly(params, Format::EVALUATION);
    // twopoly = 2;  // 这应该创建一个在所有模数上都是2的常数多项式

    // // 设置g[0] = 1
    // g(0, 0) = onepoly;

    // // 通过乘以2生成后续元素
    // DCRTPoly current = onepoly;  // 从1开始
    // for (size_t i = 1; i < k; ++i) {
    //     current = current * twopoly;  // 每次乘以2
    //     current.SetFormat(Format::EVALUATION);  // 确保格式正确
    //     g(0, i) = current;
    // }
    gadget.SetFormat(Format::EVALUATION);
        debug::printAsNormalPoly(gadget, "gadget", 10);

    // DCRTPoly g_inv = g.MultiplicativeInverse();
    //     debug::printAsNormalPoly(g_inv*g, "g_inv*g", 10);
    // // 设置格式为EVALUATION以进行后续计算



///////////////////////////////////////////

    // 计算 a2 = h*g - a'^T * T
    Matrix<DCRTPoly> a2_1 = (h * gadget.Transpose());
        debug::printAsNormalPoly(a2_1, "a2_1");
    Matrix<DCRTPoly> a2_2 = (a_prime.Transpose() * T);
        debug::printAsNormalPoly(a2_2, "a2_2");
    Matrix<DCRTPoly> a2 = (a2_1 - a2_2).Transpose();
        debug::printAsNormalPoly(a2, "a2");

        
    // 垂直拼接 a' 和 a2 得到最终的a
    // a = [a'^T| hg-a'^TT]^T
    Matrix<DCRTPoly> a(zero_alloc, a_prime.GetRows() + a2.GetRows(), a_prime.GetCols());
        // std::cout << "a dimensions: " << a.GetRows() << " x " << a.GetCols() << std::endl;
        debug::printAsNormalPoly(a, "a");
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
        // debug::printAsNormalPoly(a, "a");
            Matrix<DCRTPoly> a_part1 = Matrix<DCRTPoly>(zero_alloc, a_prime.GetRows(), a_prime.GetCols());
            a_part1.SetFormat(Format::EVALUATION);
            for (size_t i = 0; i < a_prime.GetRows(); ++i) {
                for (size_t j = 0; j < a_prime.GetCols(); ++j) {
                    a_part1(i, j) = a(i, j);
                }
            }
            debug::printAsNormalPoly(a_part1, "a_来自a_prime的部分");
            
            // 查看a矩阵中来自a2的部分（后a2.GetRows()行）
            Matrix<DCRTPoly> a_part2 = Matrix<DCRTPoly>(zero_alloc, a2.GetRows(), a2.GetCols());
            a_part2.SetFormat(Format::EVALUATION);
            for (size_t i = 0; i < a2.GetRows(); ++i) {
                for (size_t j = 0; j < a2.GetCols(); ++j) {
                    a_part2(i, j) = a(a_prime.GetRows() + i, j);
                }
            }
            debug::printAsNormalPoly(a_part2, "a_来自a2的部分");

    return std::make_pair(a, T);
}

template <class Element>
Matrix<Element> MyRLWETrapdoorUtility<Element>::SamplePre(
    const Matrix<Element>& T,
    const Matrix<Element>& a, const Element& h,
    double eta, double sigma, double alpha,
    const Element& u,
    int64_t base) {
    
    const std::shared_ptr<typename Element::Params> params = u.GetParams();
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    
    // 步骤1: Sample y ← SampleY(q, η, α, T)
        std::cout << "Begin sampling Y: " << std::endl;
    Matrix<Element> y = SampleY(params->GetModulus(), eta, alpha, T, params);
        debug::printAsNormalPoly(y, "y");
        debug::printAsNormalPoly(a, "a");
        debug::printAsNormalPoly(u, "u");

    Matrix<Element> ay = a.Transpose() * y;
        debug::printAsNormalPoly(ay, "a_times_y");
    Element a_times_y = ay(0,0);
        debug::printAsNormalPoly(a_times_y, "a_times_y");
    Element u_minus_ay = u - a_times_y;
        // debug::printAsNormalPoly(u_minus_ay, "u - a_times_y");
        debug::printAsNormalPoly(u_minus_ay + a_times_y - u, "u - a_times_y + a_times_y-u");
    // 计算h的模逆
    Element h_inv = h.MultiplicativeInverse();
        debug::printAsNormalPoly(h_inv * h, "hhinv");

    Element v = h_inv * u_minus_ay;
        debug::printAsNormalPoly(v, "v");
        // Element verification = v * h + a_times_y(0,0);
        // debug::printAsNormalPoly(h * v-u_minus_ay, "hv-u_minus_ay");
        debug::printAsNormalPoly(v * h + a_times_y - u, "v * h + a_times_y-u");

    // 步骤3: Sample z ← Sample PolyG(σ, v)
    size_t k_dim = T.GetCols();  // z的维度应与T的列数相同
    Matrix<Element> z = SamplePolyG(sigma, v, k_dim, params);
        debug::printAsNormalPoly(z, "z");
        Matrix<Element> gadget = GenGadget(k_dim, params);
        debug::printAsNormalPoly(gadget, "gadget");
        Matrix<Element> gz = gadget.Transpose() * z;
        debug::printAsNormalPoly(gz, "gz");
        debug::printAsNormalPoly(gz(0,0)-v, "gz-v");
        // Element z_sum = zero_alloc();
        
        // z_sum.SetFormat(Format::EVALUATION);
        // for (size_t i = 0; i < k_dim; ++i) {
        //     z_sum = z_sum + z(i, 0);
        // }
        // debug::printAsNormalPoly(z_sum-v, "z_sum");

    
    // 步骤4: 计算 x ← y + [T; I] · z
    size_t m_prime = T.GetRows();  // T的行数
    // size_t k_dim = T.GetCols();    // T的列数
    
    // 构建 [T; I] 矩阵
        std::cout << "Create T_I " << std::endl;
    Matrix<Element> T_I(zero_alloc, m_prime + k_dim, k_dim);
    
    // 上半部分是T
    for (size_t i = 0; i < m_prime; ++i) {
        for (size_t j = 0; j < k_dim; ++j) {
            T_I(i, j) = T(i, j);
        }
    }
        std::cout << "Fill with T " << std::endl;

    Element unity_poly = CreateUnityPoly(params);
    debug::printAsNormalPoly(unity_poly, "unity_poly");
    for (size_t i = 0; i < k_dim; ++i) {
        for (size_t j = 0; j < k_dim; ++j) {
            if (i == j) {
                T_I(m_prime + i, j) = unity_poly;
            } else {
                T_I(m_prime + i, j) = zero_alloc();
            }
        }
    }
        // a = [a'^| hg-a'^ T]^T
        // at * TI = hg
        Matrix<Element> a_times_T_I = a.Transpose() * T_I;
        debug::printAsNormalPoly(a_times_T_I, "a_times_T_I");
        // z = h-1 * (u - ay)
    Matrix<Element> T_I_times_z = T_I * z;
        debug::printAsNormalPoly(T_I, "TI");
        debug::printAsNormalPoly(T_I_times_z, "T_I_times_z");
        
        // std::cout << "T_I·z dimensions: " << T_I_times_z.GetRows() << " x " << T_I_times_z.GetCols() << std::endl;
    Matrix<Element> x = y + T_I_times_z;
        debug::printAsNormalPoly(x, "x");
    
    return x;
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
    // 创建一个零多项式
    auto alloc = Element::Allocator(params, COEFFICIENT);
    Element one_poly = alloc();
    
    // 对于DCRTPoly，使用不同的方法
    if constexpr (std::is_same<Element, DCRTPoly>::value) {

        DCRTPoly one_dcrt(params, Format::EVALUATION);
        one_dcrt = 1;  // 设置为常数
        one_poly = one_dcrt; // true表示初始化为1
    } else {
        // 对于Poly和NativePoly，使用原来的方法
        typename Element::Vector coefficients = one_poly.GetValues();
        if (coefficients.GetLength() > 0) {
            coefficients[0] = 1;
        }
        one_poly.SetValues(coefficients, COEFFICIENT);
    }
    
    one_poly.SetFormat(EVALUATION);
    return one_poly;
}

template <class Element>
Matrix<Element> MyRLWETrapdoorUtility<Element>::SamplePolyG(
    double sigma, const Element& v, 
    size_t k, std::shared_ptr<typename Element::Params> params) {
    
    Matrix<Element> gadget = GenGadget(k, params);
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    
    // 计算 α = √5 * σ
    double sampling_sigma = std::sqrt(5.0) * sigma;
    
    // 创建高斯分配器
    auto gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(
        params, COEFFICIENT, sampling_sigma);
    
    // 采样z，维度为k x 1
    Matrix<Element> z(zero_alloc, k, 1, gaussian_alloc);
    z.SetFormat(Format::EVALUATION);
    
    // for (size_t i = 0; i < k-1; ++i) {
    //     Element poly = gaussian_alloc();
    //     poly.SetFormat(Format::EVALUATION);
    //     z(i, 0) = poly;
    // }
    Element weighted_sum = zero_alloc();
    weighted_sum.SetFormat(Format::EVALUATION);
    for (size_t i = 0; i < k-1; ++i) {
        weighted_sum = weighted_sum + gadget(i, 0) * z(i, 0);
    }

    Element diff = v - weighted_sum;
    diff.SetFormat(Format::EVALUATION);
    
    // 计算gadget向量的最后一个元素 g[k-1] 的模逆
    Element g_last = gadget(k-1, 0);
    Element g_last_inv = g_last.MultiplicativeInverse();
    Element z_last = diff * g_last_inv;
    z_last.SetFormat(Format::EVALUATION);
    z(k-1, 0) = z_last;
    
    // for (size_t i = 0; i < k-1; ++i) {
    //     z_sum = z_sum + z(i, 0);
    // }
    // Element diff = v - z_sum;
    // diff.SetFormat(Format::EVALUATION);
    // z(k-1, 0) = diff;
    // z_sum = z_sum + diff;
    // debug::printAsNormalPoly(z_sum, "z_sum_in");

    return z;
}

template <class Element>
Matrix<Element> MyRLWETrapdoorUtility<Element>::SampleY(
    const typename Element::Integer& q, double eta, double alpha,
    const Matrix<Element>& T, 
    std::shared_ptr<typename Element::Params> params) {
    
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    
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
Matrix<Element> MyRLWETrapdoorUtility<Element>::GenGadget(const size_t k, std::shared_ptr<typename Element::Params> params) 
{
    // 分配零矩阵
    auto zero_alloc = Element::Allocator(params, EVALUATION);
    Matrix<Element> g = Matrix<Element>(zero_alloc, k, 1);
    
    // 创建常数多项式
    Element onepoly(params, Format::EVALUATION);
    onepoly = 1;  // 设置为常数1
    
    Element twopoly(params, Format::EVALUATION);
    twopoly = 2;  // 设置为常数2
    
    // 设置g[0] = 1
    g(0, 0) = onepoly;
    
    // 通过乘以2生成后续元素
    Element current = onepoly;  // 从1开始
    for (size_t i = 1; i < k; ++i) {
        current = current * twopoly;  // 每次乘以2
        current.SetFormat(Format::EVALUATION);  // 确保格式正确
        g(i, 0) = current;
    }
    
    return g;
}

template <class Element>
Matrix<Element> MyRLWETrapdoorUtility<Element>::GenSecretKey(size_t m_minus_k, double sigma, std::shared_ptr<typename Element::Params> params)
{
    auto zero_alloc = Element::Allocator(params, Format::EVALUATION);
    auto gaussian_alloc = Element::MakeDiscreteGaussianCoefficientAllocator(params, Format::COEFFICIENT, sigma);
    
    Matrix<Element> s(zero_alloc, m_minus_k, 1, gaussian_alloc);
    s.SetFormat(Format::EVALUATION);

    return s;
}




// 显式实例化常用类型
template class MyRLWETrapdoorUtility<Poly>;
template class MyRLWETrapdoorUtility<NativePoly>;
template class MyRLWETrapdoorUtility<DCRTPoly>;
} // namespace myrlwe

