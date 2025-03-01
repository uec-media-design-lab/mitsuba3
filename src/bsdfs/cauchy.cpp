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
class cauchy final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    using Index = dr::uint32_array_t<Float>;
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    cauchy(const Properties &props) : Base(props)
    {
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        m_divergence = props.texture<Texture>("divergence", 1.f);
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        m_a = props.texture<Texture>("a", 1.f);

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

    void traverse(TraversalCallback *callback) override
    {
        callback->put_object("divergence", m_divergence.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);

        isTraversed = true;
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
    {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float sin_theta_i = Frame3f::sin_theta(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float sin_phi_i = Frame3f::sin_phi(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};


        // ----- 出射光のサンプリング -----
        // theta, phiのサンプリング
        Float del_phi = 2.0f*dr::Pi<Float> * sample2.x(); // 角度の差分にする
        Float del_theta = m_a->eval(si, active).x() * dr::tan(dr::Pi<Float> * (sample2.y() - 0.5f));
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_theta);
        Vector3f norm = si.n;
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;

        bs.sampled_component = 0;
        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        bs.eta = 1.f;
        bs.pdf = 1.f;

        UnpolarizedSpectrum cauchyvalue = 1.0f * dr::InvPi<Float> / (1 + dr::pow(del_theta, 2));
        cauchyvalue = dr::select(
            dr::isfinite(cauchyvalue),
            cauchyvalue,
            ScalarFloat(1.0)
        );

        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * cauchyvalue;

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

        UnpolarizedSpectrum cauchyvalue = 1.0f * dr::InvPi<Float> / (1 + dr::pow(delta_theta, 2));
        cauchyvalue = dr::select(
            dr::isfinite(cauchyvalue),
            cauchyvalue,
            ScalarFloat(1.0)
        );
        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * cauchyvalue;

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

        UnpolarizedSpectrum cauchyvalue = 1.0f * dr::InvPi<Float> / (1 + dr::pow(delta_theta, 2));
        cauchyvalue = dr::select(
            dr::isfinite(cauchyvalue),
            cauchyvalue,
            ScalarFloat(1.0)
        );

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, cauchyvalue[0], 0.f);
    }

    std::pair<Spectrum, Float> eval_pdf(
        const BSDFContext &ctx,
        const SurfaceInteraction3f &si,
        // SurfaceInteraction3f &si,
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

            UnpolarizedSpectrum cauchyvalue = 1.0f * dr::InvPi<Float> / (1 + dr::pow(delta_theta, 2));
            cauchyvalue = dr::select(
                dr::isfinite(cauchyvalue),
                cauchyvalue,
                ScalarFloat(1.0)
            );
            UnpolarizedSpectrum value = m_reflectance->eval(si, active) * cauchyvalue;
            Float pdf = dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, cauchyvalue[0], 0.f);

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
        }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "cauchy[" << std::endl
            << " divergence = " << string::indent(m_divergence) << std::endl
            << " reflectance = " << string::indent(m_reflectance) << std::endl
            << "a = " << string::indent(m_a) << std::endl
            << "]";
        return oss.str();
    }


    MI_DECLARE_CLASS()

// フィールド
private:
    ref<Texture> m_divergence;
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    bool isTraversed;
};

MI_IMPLEMENT_CLASS_VARIANT(cauchy, BSDF)
MI_EXPORT_PLUGIN(cauchy, "cauchy")
NAMESPACE_END(mitsuba)