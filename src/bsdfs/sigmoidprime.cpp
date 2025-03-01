#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/texture.h>

#include <drjit/dynamic.h>
#include <mitsuba/core/distr_2d.h>
#include <mitsuba/core/distr_1d.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <random>
#include <iterator>
#include <string>
// #include "nlohmann/json.hpp"
using namespace std;
// using json = nlohmann::json;


NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class sigmoidprime2 final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    using Index = dr::uint32_array_t<Float>;
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    sigmoidprime2(const Properties &props) : Base(props)
    {
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        m_divergence = props.texture<Texture>("divergence", 1.f);
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        m_a = props.texture<Texture>("a", 1.f);

        // M = 720; N = 180;
        // M = 1440; N = 360;
        // M = 2880; N = 720;        
        M = 1; N = 30000;
        // M = 4320; N = 1080;
        // N = 50000;

        // offset = dr::Pi<Float>/15555.0f;
        // offset = dr::Pi<Float>/10000.0f;
        // offset = dr::Pi<Float>/20000.0f;
        // offset = dr::Pi<Float>/2000.0f;
        offset = 0.f;

        isTraversed = false;
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

    

    void make_sigmoidprime1(float a, int expo)
    {
        ScalarFloat data[N];
        int i = 0, idx = 0;

        for (int i = 0; i < N; ++i)
        {
            // ScalarFloat arg = a * (i - N / 2.0f);
            // ScalarFloat arg = dr::Pi<Float> * i / N;
            ScalarFloat arg = a * (dr::Pi<Float>*i / (2.0f*N));
            // ScalarFloat arg = a * (20.0f*i / N);
            ScalarFloat sigmoidvalue = 1.0f / (1 + dr::exp(-arg));
            ScalarFloat sincvalue = sigmoidvalue*(1-sigmoidvalue);
            // Mask mask = arg == 0.0;
            // argが0の時はゼロ除算を防ぐため
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
        }
        // std::cout << "[" << std::endl;
        // float sum = dr::Pi<Float>;
        for (int i = 0; i < N; ++i)
        {
            data[i] /= sum;
            // data[i] *= 0.4;
        }

        sum_rawPDF = sum;
        m_data = &data[0]; // LUTの先頭ポインタ
        struct DiscreteDistribution<Float> dd(m_data, N);
        m_pdfdata = dd;
    }

    void make_sigmoidprime1_multi(float a, int expo)
    {
        ScalarFloat data[N];
        int i = 0, idx = 0;

        for (int i = 0; i < N; ++i)
        {
            ScalarFloat sincvalue = ScalarFloat(0.0f);
            for (int n = -100; n <= 100; ++n) {
                ScalarFloat arg = a * ((dr::Pi<Float>*i / (2.0f*N)) - n*0.08);
                sincvalue += dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo);
                // argが0の時はゼロ除算を防ぐため
                sincvalue = dr::select(
                    dr::isfinite(sincvalue),
                    dr::abs(sincvalue),
                    ScalarFloat(1.0)
                );
            }
            
            data[i] = sincvalue;
        }

        // 正規化（ごり押し）
        float sum = 0.0f;
        for (int i = 0; i < N; ++i)
        {
            sum += data[i];
        }
        // float sum = dr::Pi<Float>;
        for (int i = 0; i < N; ++i)
        {
            data[i] /= sum;
            // data[i] *= 0.4;
        }

        sum_rawPDF = sum;
        m_data = &data[0]; // LUTの先頭ポインタ
        struct DiscreteDistribution<Float> dd(m_data, N);
        m_pdfdata = dd;
    }

    void traverse(TraversalCallback *callback) override
    {
        callback->put_object("divergence", m_divergence.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);

        isTraversed = true;

        // mi.traverse(scene)が呼び出されたときに新しいaのLUTを作る
        // M = 1440; N = 360;
        // make_sigmoidprime2(m_a.get()->max(), 2);

        make_sigmoidprime1(m_a.get()->max(), 2);
        // make_sigmoidprime1_multi(m_a.get()->max(), 2);

        Float a = m_a.get()->max();
        // offset2 = (a+2)*dr::Pi<Float>/(a*(a+1));
        // offset2 = dr::Pi<Float> / a;
        offset2 = 0.f;
    }

    float generate_random() const
    {
        random_device seed_gen;
        static std::mt19937 generator(seed_gen());
        static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
        return distribution(generator);
    }


    std::pair<BSDFSample3f, Spectrum> sample(
        const BSDFContext &ctx,
        const SurfaceInteraction3f &si,
        // SurfaceInteraction3f &si,
        Float sample1,
        const Point2f &sample2,
        Mask active
        ) const override
        // ) const
    {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float sin_theta_i = Frame3f::sin_theta(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float sin_phi_i = Frame3f::sin_phi(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};

        // Point2u point2;
        // Point2f point2;
        Float point;
        Float pdfvalue;
        Point2f sampled;
        // std::tie(point2, pdfvalue, sampled) = m_pdfdata2.sample(sample2, active);
        // point = m_pdfdata.sample(sample1, active);
        point = m_pdfdata.sample(sample2.y());

        // ----- 出射光のサンプリング -----
        // theta, phiのサンプリング
        // index_theta = point;

        Float del_phi = dr::TwoPi<Float> * sample2.x(); // 角度の差分にする
        printf("phi = %f\n", (del_phi*180.f)/dr::Pi<Float>);
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N) + dr::Pi<Float>/15555.0f;
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N) + dr::Pi<Float>/(2.0*m_a->eval(si, active).x());
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N) + dr::Pi<Float>*0.0276583f;
        Float a = m_a->eval(si, active).x();
        Float del_theta = point*dr::Pi<Float>/(2.0f*N) + offset + offset2;
        // Float tan_theta = point*dr::Pi<Float>/(2.0f*N);
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N);
        // Float del_theta = dr::atan(tan_theta);
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N);
        // printf("%f\n", point);
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_theta);
        // Vector3f delvec = Vector3f(0.f, 0.f, 1.f);
        // Vector3f norm = si.n;
        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        // printf("%f, %f, %f\n", si.wi.x(), si.wi.y(), si.wi.z());
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;

        bs.sampled_component = 0;
        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        bs.eta = 1.f;
        bs.pdf = 1.f;

        // SurfaceInteraction3f& si_nonconst = const_cast<SurfaceInteraction3f&>(si);
        // si_nonconst.p += Point3f(0.065*cos_del_phi, 0.065*sin_del_phi, 0.0);

        UnpolarizedSpectrum arg = m_a->eval(si, active)*del_theta;
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask mask = arg == 0;
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue /= sum_rawPDF;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active);

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

        Float a = m_a->eval(si, active).x();
        UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset - offset2);
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask mask = arg == 0;
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue /= sum_rawPDF;
        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;

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

        Float a = m_a->eval(si, active).x();
        // Float offset2 = 2*dr::Pi<Float>/(a*(a+1));
        UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset - offset2);
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
        // UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask mask = arg == 0;
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue /= sum_rawPDF;

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue.x(), 0.f);
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

            Float a = m_a->eval(si, active).x();
            // Float offset = dr::Pi<Float>/15555.0f;
            // Float offset2 = 2*dr::Pi<Float>/(a*(a+1));
            UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset - offset2);
            UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
            // UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
            // UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
            // Mask mask = arg == 0;
            sincvalue = dr::select(
                dr::isfinite(sincvalue),
                sincvalue,
                ScalarFloat(1.0)
            );
            // sincvalue /= sum_rawPDF;
            UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue;
            Float pdf = dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue.x(), 0.f);

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
        }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "sigmoidprime2[" << std::endl
            << " divergence = " << string::indent(m_divergence) << std::endl
            << " reflectance = " << string::indent(m_reflectance) << std::endl
            << "a = " << string::indent(m_a) << std::endl
            << "]";
        return oss.str();
    }

    // std::string to_string() const override {
    //     std::ostringstream oss;
    //     oss << "a = " << m_a.get() << " .json" << std::endl;
    //     cout << "a = " << m_a.get() << std::endl;
    //     return oss.str();
    // }

    std::string make_LUTPath() const {
        std::ostringstream oss;
        std::ostringstream out;
        int precision = 2;
        out << std::fixed << std::setprecision(precision) << m_a.get()->max();
        std::string a_str = out.str();
        oss << LUTDir << "pdf_sinc2, a = " << a_str << ".json";
        return oss.str();
    }


    MI_DECLARE_CLASS()

// フィールド
private:
    uint32_t M; // LUTのphiサイズ
    uint32_t N; // LUTのthetaサイズ
    float sum_rawPDF;
    Point2u point2;
    ref<Texture> m_divergence;
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    vector<vector<float>> m_beforedata;
    ScalarFloat *m_data;
    DiscreteDistribution<Float> m_pdfdata;
    DiscreteDistribution2D<Float, 2> m_pdfdata2;
    // Hierarchical2D<Float, 2> m_pdfdata;
    ScalarVector2u m_size;
    std::string LUTFilename;
    std::string LUTDir;
    Float index_phi, index_theta;
    bool isTraversed;
    Float offset;
    Float offset2;
};

MI_IMPLEMENT_CLASS_VARIANT(sigmoidprime2, BSDF)
MI_EXPORT_PLUGIN(sigmoidprime2, "sigmoidprime2")
NAMESPACE_END(mitsuba)