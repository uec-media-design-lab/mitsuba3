#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/texture.h>

#include <drjit/dynamic.h>
#include <mitsuba/core/distr_2d.h>
#include <mitsuba/core/distr_1d.h>

#include <mitsuba/render/sampler.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <random>
#include <iterator>
#include <string>
#include "nlohmann/json.hpp"
using namespace std;
using json = nlohmann::json;


NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class sinc_shift3 final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    using Index = dr::uint32_array_t<Float>;
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    sinc_shift3(const Properties &props) : Base(props)
    {
        m_rrrw = m_rrrh = 175.0;
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        m_a = props.texture<Texture>("a", 1.f);
        m_angle = props.texture<Texture>("angle", 0.f);
        m_cornersize = props.texture<Texture>("cornersize", 0.f);
        m_shiftoffset = props.texture<Texture>("shiftoffset", 0.f);

        M = 1; N = 30000;

        // offset = dr::Pi<Float>/15555.0f;
        // offset = 0.0001555555555554644*dr::Pi<Float>;
        // offset = 0.0003*dr::Pi<Float>;
        // offset = 18.0/(3600.0*180.0)*dr::Pi<Float>;
        offset = 0.0;

        // make_rayshift_array();

        m_isEstimation = 1;
    }

    // class - rotations
    Float mag(Normal3f n) const {
        return dr::sqrt(n.x()*n.x() + n.y()*n.y() + n.z()*n.z());
    }

    Vector3f rotate (const Vector3f vec, const Vector3f normal, const Vector3f to) const
    {
        Float cos_theta_i = Frame3f::cos_theta(to);
        Float sin_theta_i = Frame3f::sin_theta(to);
        Float cos_phi_i = Frame3f::cos_phi(to);
        Float sin_phi_i = Frame3f::sin_phi(to);

        Vector3f axis = dr::normalize(dr::cross(normal, to));

        // ロドリゲスの回転公式
        Vector3f fromNorm = dr::normalize(normal);
        Vector3f toNorm = dr::normalize(to);
        Float s = mag(dr::cross(fromNorm, toNorm));
        Float c = dr::dot(fromNorm, toNorm);
        Vector3f result = c*vec + dr::dot(axis, vec)*(1-c)*axis + dr::cross(axis, vec)*s;
        result = dr::normalize(result);

        return result;
    }

    void make_sincarray1(float a, int expo, int angle)
    {
        ScalarFloat data[N];
        int i = 0, idx = 0;

        for (int i = 0; i < N; ++i)
        {
            ScalarFloat arg = (a*300/dr::cos(offset)) * dr::tan((dr::Pi<Float>*i/(2.0f*N) - offset));
            ScalarFloat sincvalue = dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo);

            // if (i >= 28000) sincvalue = 0.f;
            // printf("%d = %f ", i, sincvalue);
            sincvalue = dr::select(
                dr::isfinite(sincvalue),
                dr::abs(sincvalue),
                ScalarFloat(1.0)
            );
            data[i] = sincvalue;
        }
        // 正規化（ごり押し）
        float sum = 0.0f;
        for (int i = 0; i < N; ++i)
        {
            sum += data[i];
            // printf("data = %f\n", data[i]);
        }
        for (int i = 0; i < N; ++i)
        {
            data[i] /= sum;
        }

        sum_rawPDF = sum;
        m_data = &data[0]; // LUTの先頭ポインタ
        struct DiscreteDistribution<Float> dd(m_data, N);
        m_pdfdata = dd;
    }

    // TODO 角度ごとのLUTを保存するように書き替える
    void make_rayshift_array() {
        vector<vector<float>> m_rayshift_list;
        std::string filename = "/home/sugawara.ryo/makeLUT/2025.03.17diff/LUTs/rayshift_LUT.json";
        ifstream ifs(filename.c_str());
        if (ifs.good())
        {
            json m_json;
            ifs >> m_json;

            for (const auto& json : m_json) {
                vector<float> d = json.get<vector<float>>();
                m_rayshift_list.push_back(d);
            }
        }
        else
        {
            cout << "ファイルの読み込みに失敗しました" << endl;
        }

        // printf("doing\n");
        angle_index = angle / 5;
        int num_i = m_rayshift_list.size();
        // int num_j = m_rayshift_list[0].size();
        int num_j = m_rayshift_list[angle_index].size();
        // printf("angle_index = %d\n", angle_index);
        // printf("data[idx] = %f\n", m_rayshift_list[angle_index][0]);
        ScalarFloat data[num_j];
        int i = 0, idx = 0;
        for (int j = 0; j < num_j; ++j) {
            data[idx++] = m_rayshift_list[angle_index][j];
        }
        // printf("\n");
        printf("rayshift_array[0] = %f, angle = %d, angle_index = %d\n", data[0], angle, angle_index);

        ScalarVector2u size(num_j, num_i);
        struct DiscreteDistribution<Float> dd(&data[0], num_j);
        m_rayshiftLUT = dd;
        rayshift_N = num_j;
        printf("rayshift_N = %d\n", num_j);
    }

    void traverse(TraversalCallback *callback) override
    {
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("angle", m_angle.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("cornersize", m_cornersize.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("shiftoffset", m_shiftoffset.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        // callback->put_parameter("isEstimation", m_isEstimation);

        // printf("isEstimation = %d\n", m_isEstimation);

        // float a = m_a.get()->max();
        angle = static_cast<int>(m_angle.get()->max());
        // printf("angle = %d\n", angle);

        angle_index = angle / 5.0;
        // float a, reflectance;
        if (m_isEstimation) {
            a = m_a.get()->max(); reflectance = m_reflectance.get()->max();
        } else {
            a = alist[angle_index]; reflectance = reflist[angle_index];
        }
        
        // make_sincarray1(m_a.get()->max(), 2, angle);
        make_sincarray1(a, 2, angle);
        make_rayshift_array();
    }

    float generate_random() const
    {
        random_device seed_gen;
        static std::mt19937 generator(seed_gen());
        static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
        return distribution(generator);
    }


    mitsuba::PCG32<UInt32> setRandomGenerator(const Float seed) const {
        mitsuba::PCG32<UInt32> rng(1, PCG32_DEFAULT_STATE, seed);
        rng.state = seed*1000000;
        rng.template next_float<Float>();
        return rng;
    }
    Float rand(mitsuba::PCG32<UInt32> &generator) const {
        return generator.template next_float<Float>();
    }
    Point2f rand2(mitsuba::PCG32<UInt32> &generator) const {
        Point2f sample2;
        sample2.x() = rand(generator);
        sample2.y() = rand(generator);
        return sample2;
    }


    std::pair<BSDFSample3f, Spectrum> sample(
        const BSDFContext &ctx,
        const SurfaceInteraction3f &si,
        Float sample1,
        const Point2f &sample2,
        Mask active
        ) const override
    {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float sin_theta_i = Frame3f::sin_theta(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float sin_phi_i = Frame3f::sin_phi(si.wi);
        Float theta_i = dr::acos(cos_theta_i);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};

        auto rng = setRandomGenerator(sample1*10000000);

        // BSSRDF
        Point3f p = si.p; Point2f r2f = rand2(rng);
        Float sampled_p2f; Float pdfvalue; Point2f sampled;
        Float r1 = rand(rng); Float r2 = rand(rng);
        // sampled_p2f = m_rayshiftLUT.sample(r1);
        sampled_p2f = m_rayshiftLUT.sample(sample1);

        float cornersize = m_cornersize.get()->max();
        float shiftoffset = m_shiftoffset.get()->max();
        Float shift = (cornersize*sampled_p2f*1.50f)/rayshift_N; // TODO あとでちゃんとここを作る
        Float r3 = rand(rng);
        Float theta = 2*r3*dr::Pi<Float>;

        bs.p = Point3f(shift*dr::cos(theta), shift*dr::sin(theta), shiftoffset);

        // ----- 出射光のサンプリング -----
        Float point;
        point = m_pdfdata.sample(sample2.y());
        // printf("point = %f\n", point);
        point = dr::select(point >= 0.f, point, 0.f); // point < 0のとき0にする
        Float random2 = rand(rng), random3 = rand(rng);
        Float del_phi = 2.0*dr::Pi<Float> * random2; // 角度の差分にする

        Float random1 = rand(rng);
        Float offset2 = 0.0001555555555554644*dr::Pi<Float>;
        Float del_theta = dr::Pi<Float>*point/(2.0*N) + offset2;
        // Float del_theta = dr::Pi<Float>*point/(2.0*N);

        Float del_phi_offset = (2.0f*dr::Pi<Float>)*random3; // 角度の差分にする
        Float cos_offset_phi = dr::cos(del_phi_offset), sin_offset_phi = dr::sin(del_phi_offset);
        Float cos_offset_theta = dr::cos(offset), sin_offset_theta = dr::sin(offset);
        Vector3f offsetvec = Vector3f(sin_offset_theta*cos_offset_phi, sin_offset_theta*sin_offset_phi, cos_offset_theta);        
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_theta);

        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;

        // Vector3f wo = rotate(delvec, norm, offsetvec);
        // Vector3f wo_true = rotate(wo, norm, si.wi);
        // bs.wo = wo_true;

        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        bs.eta = 1.f;

        // Float sincvalue = m_pdfdata.eval_pmf_floatindex(point, active)*sum_rawPDF;
        // UnpolarizedSpectrum value = reflectance*sincvalue;
        // bs.pdf = sincvalue;
        
        // float a, reflectance;
        // if (m_isEstimation) {a = m_a.get()->max(); reflectance = m_reflectance.get()->max();}
        // else {a = alist[angle_index]; reflectance = reflist[angle_index];}
        // float reflectance = m_reflectance.get()->max();

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);
        bs.pdf = 1.f;
        // printf("value = %f\n", value.x());

        return {bs, depolarizer<Spectrum>(value) & (active && bs.pdf > 0.f)};
    }

    Spectrum eval(
        const BSDFContext &ctx,
        const SurfaceInteraction3f &si,
        const Vector3f &wo,
        Mask active
    ) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);
        Float sin_theta_i = Frame3f::sin_theta(si.wi),
              sin_theta_o = Frame3f::sin_theta(wo);

        Float cos_delta_theta = cos_theta_o*cos_theta_i + sin_theta_o*sin_theta_i;
        Float delta_theta = acos(cos_delta_theta);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        // Float arg = (m_a->eval(si, active).x()*300/dr::cos(offset)) * dr::tan(delta_theta - offset);

        float a, reflectance;
        // if (m_isEstimation) {a = m_a.get()->max(); reflectance = m_reflectance.get()->max();}
        // else {a = alist[angle_index]; reflectance = reflist[angle_index];}
        a = m_a.get()->max(); reflectance = m_reflectance.get()->max();
        
        Float arg = (a*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
        Float sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            ScalarFloat(1.0),
            sincvalue
        );
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        // UnpolarizedSpectrum value = reflectance * sincvalue;

        return depolarizer<Spectrum>(value) & active;
    }

    Float pdf(
        const BSDFContext &ctx,
        const SurfaceInteraction3f &si,
        const Vector3f &wo,
        Mask active
    ) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;
        
        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);
        Float sin_theta_i = Frame3f::sin_theta(si.wi),
              sin_theta_o = Frame3f::sin_theta(wo);

        Float cos_delta_theta = cos_theta_o*cos_theta_i + sin_theta_o*sin_theta_i;
        Float delta_theta = acos(cos_delta_theta);

        float a, reflectance;
        // if (m_isEstimation) {a = m_a.get()->max(); reflectance = m_reflectance.get()->max();}
        // else {a = alist[angle_index]; reflectance = reflist[angle_index];}
        a = m_a.get()->max(); reflectance = m_reflectance.get()->max();

        Float arg = (a*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
        Float sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            ScalarFloat(1.0),
            sincvalue
        );
        return sincvalue;
    }

    std::pair<Spectrum, Float> eval_pdf(
        const BSDFContext &ctx,
        const SurfaceInteraction3f &si,
        const Vector3f &wo, Mask active
        ) const override {
            MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

            if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
                return {0.f, 0.f};
            
            Float cos_theta_i = Frame3f::cos_theta(si.wi),
                  cos_theta_o = Frame3f::cos_theta(wo);
            Float sin_theta_i = Frame3f::sin_theta(si.wi),
                  sin_theta_o = Frame3f::sin_theta(wo);

            Float cos_delta_theta = cos_theta_o*cos_theta_i + sin_theta_o*sin_theta_i;
            Float delta_theta = acos(cos_delta_theta);

            active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

            // Float arg_Float = (m_a->eval(si, active).x()*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
            
            float a, reflectance;
            // if (m_isEstimation) {a = m_a.get()->max(); reflectance = m_reflectance.get()->max();}
            // else {a = alist[angle_index]; reflectance = reflist[angle_index];}
            a = m_a.get()->max(); reflectance = m_reflectance.get()->max();

            Float arg_Float = (a*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
            Float sincvalue_pdf = dr::pow(dr::sin(arg_Float), 2) / dr::pow(arg_Float, 2);
            sincvalue_pdf = dr::select(
                dr::isnan(sincvalue_pdf),
                ScalarFloat(1.0),
                sincvalue_pdf
            );
            UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue_pdf;
            // UnpolarizedSpectrum value = m_reflectance->eval(si, active);
            // UnpolarizedSpectrum value = reflectance;
            Float pdf = 1.f;

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
        }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "sinc_shift3[" << std::endl
            << " reflectance = " << string::indent(m_reflectance) << std::endl
            << "a = " << string::indent(m_a) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

// フィールド
private:
    uint32_t M; // LUTのphiサイズ
    uint32_t N; // LUTのthetaサイズ
    uint32_t rayshift_N; // rayshift_LUTのphiサイズ
    float sum_rawPDF;
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    float a, reflectance;
    int angle, angle_index;
    ScalarFloat *m_data;
    DiscreteDistribution<Float> m_pdfdata;
    DiscreteDistribution<Float> m_rayshiftLUT;
    ScalarFloat offset;
    ref<Texture> m_angle;
    ref<Texture> m_cornersize;
    ref<Texture> m_shiftoffset;
    Float m_rrrw;
    Float m_rrrh;
    bool m_isEstimation;
    float alist[10] = {9.93020, 7.88426, 7.84838, 7.15263, 7.38482, 5.10610, 3.51115, 2.45607, 1.68079, 1.17931};
    float reflist[10] = {0.33505, 0.296091, 0.28343, 0.26988, 0.26051, 0.21822, 0.16955, 0.15701, 0.15030, 0.10973};
};

MI_IMPLEMENT_CLASS_VARIANT(sinc_shift3, BSDF)
MI_EXPORT_PLUGIN(sinc_shift3, "sinc_shift3")
NAMESPACE_END(mitsuba)