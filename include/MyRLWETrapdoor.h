// MyRLWETrapdoor.h
#pragma once

#include "openfhe.h"

using namespace lbcrypto;

namespace myrlwe {

// 自定义陷门对结构
template <class Element>
class MyRLWETrapdoorPair {
public:
    Matrix<Element> T;  // 陷门矩阵
    
    MyRLWETrapdoorPair() : T(Matrix<Element>([]() { return Element(); }, 0, 0)) {}
    
    MyRLWETrapdoorPair(const Matrix<Element>& trapdoor) : T(trapdoor) {}
    
    template <class Archive>
    void serialize(Archive& ar, std::uint32_t const version) {
        ar(CEREAL_NVP(T));
    }
};

template <class Element>
class MyRLWETrapdoorUtility {
public:
    using ParmType = typename Element::Params;
    using DggType = typename Element::DggType;
    using IntType = typename Element::Integer;
    
    /**
     * 你的Trapdoor生成算法
     */
    static std::pair<Matrix<Element>, Matrix<Element>> TrapdoorGen(
        std::shared_ptr<ParmType> params, const Matrix<Element>& a_prime, 
        const Element& h, double sigma, int64_t base = 2);

    /**
     * 你的SamplePre算法
     */
    static Matrix<Element> SamplePre(const Matrix<Element>& T,
                                    const Matrix<Element>& a, const Element& h,
                                    double eta, double sigma, double alpha,
                                    const Element& u,
                                    int64_t base = 2);
    static Element CreateUnityPoly(std::shared_ptr<ParmType> params);

    static Matrix<Element> GenSecretKey(size_t m_minus_k, double sigma, std::shared_ptr<ParmType> params);
private:
    /**
     * SamplePolyG算法
     */
    static Matrix<Element> SamplePolyG(double sigma, const Element& v,
                                      size_t k,
                                      std::shared_ptr<ParmType> params);

    /**
     * SampleY算法
     */
    static Matrix<Element> SampleY(const IntType& q, double eta, double alpha,
                                  const Matrix<Element>& T,
                                std::shared_ptr<ParmType> params);

    /**
     * 估算矩阵的谱范数
     */


    static double EstimateSpectralNorm(const Matrix<Element>& T);
    static Matrix<Element> GenGadget(const size_t k, std::shared_ptr<ParmType> params);

};

} // namespace myrlwe