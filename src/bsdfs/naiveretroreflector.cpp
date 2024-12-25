#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

#include <mitsuba/core/ray.h>
#include <drjit/transform.h>
// #include <mitsuba/render/microfacet.h>
// #include <mitsuba/core/transform.h>
// #include <iostream>
// #include <fstream>
#include <random>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _bsdf-diffuse:

Smooth diffuse material (:monosp:`diffuse`)
-------------------------------------------

.. pluginparameters::

 * - reflectance
   - |spectrum| or |texture|
   - Specifies the diffuse albedo of the material (Default: 0.5)
   - |exposed|, |differentiable|

The smooth diffuse material (also referred to as *Lambertian*)
represents an ideally diffuse material with a user-specified amount of
reflectance. Any received illumination is scattered so that the surface
looks the same independently of the direction of observation.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_plain.jpg
   :caption: Homogeneous reflectance
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_textured.jpg
   :caption: Textured reflectance
.. subfigend::
   :label: fig-diffuse

Apart from a homogeneous reflectance value, the plugin can also accept
a nested or referenced texture map to be used as the source of reflectance
information, which is then mapped onto the shape based on its UV
parameterization. When no parameters are specified, the model uses the default
of 50% reflectance.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the :ref:`twosided <bsdf-twosided>` BRDF adapter plugin.
The following XML snippet describes a diffuse material,
whose reflectance is specified as an sRGB color:

.. tabs::
    .. code-tab:: xml
        :name: diffuse-srgb

        <bsdf type="diffuse">
            <rgb name="reflectance" value="0.2, 0.25, 0.7"/>
        </bsdf>

    .. code-tab:: python

        'type': 'diffuse',
        'reflectance': {
            'type': 'rgb',
            'value': [0.2, 0.25, 0.7]
        }

Alternatively, the reflectance can be textured:

.. tabs::
    .. code-tab:: xml
        :name: diffuse-texture

        <bsdf type="diffuse">
            <texture type="bitmap" name="reflectance">
                <string name="filename" value="wood.jpg"/>
            </texture>
        </bsdf>

    .. code-tab:: python

        'type': 'diffuse',
        'reflectance': {
            'type': 'bitmap',
            'filename': 'wood.jpg'
        }
*/
template <typename Float, typename Spectrum>
class naiveRetroReflector final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    naiveRetroReflector(const Properties &props) : Base(props) {
        m_reflectance = props.texture<Texture>("reflectance", .5f);
        m_alpha = props.texture<Texture>("alpha", .00005f);
        m_flags = BSDFFlags::DeltaReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);

        std::cout << "naiveRetroReflector Initialized." << std::endl;
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("reflectance", m_reflectance.get(), +ParamFlags::Differentiable);
        callback->put_object("alpha", m_alpha.get(), +ParamFlags::Differentiable);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* sample1 */,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        // std::cout << "naiveRetroReflector Sample" << std::endl;

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float theta_i_deg = dr::acos(cos_theta_i) * 180.0 * dr::InvPi<Float>;
        Float theta_i = dr::acos(cos_theta_i);
        Float phi_i_deg = dr::acos(cos_phi_i) * 180.0 * dr::InvPi<Float>;
        Float phi_i = dr::acos(cos_phi_i);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::DeltaReflection)))
            return { bs, 0.f };

        // std::cout << "cos_theta_i = " << (dr::acos(cos_theta_i)*180.0*dr::InvPi<Float>) << std::endl;
        // std::cout << "cos_phi_i = " << (dr::acos(cos_phi_i)*180.0*dr::InvPi<Float>) << std::endl;
        
        // std::normal_distribution<> dist(theta_i_deg, 0.0005);
        // float theta_o_deg = dist(sample2.y());
        // std::cout << "theta_o_deg = " << theta_o_deg << std::endl;

        // とりあえずGGXに従う乱数を生成する
        // 中心theta_iであり、標準偏差は探索するパラメータ
        // MicrofacetDistribution microfacet(MicrofacetType::GGX, m_alpha->eval_1(si, active));
        // Vector3f m;
        // std::tie(m, bs.pdf) = microfacet.sample(si.wi, sample2);

        // bs.wo = Transform.rotate(si.n, phi_i)*Transform.rotate(si.n, theta_i)*si.n;
        // static constexpr size_t Size = Point_::Size;
        // struct Transform transform = {};
        using Value = Float;
        using Matrix4f = dr::Matrix<Value, 4>;
        // bs.wo = dr::matmul(dr::matmul(dr::rotate<Matrix4f>(Vector3f(0.f, 1.f, 0.f), theta_i), dr::rotate<Matrix4f>(si.n, phi_i)), si.n);
        bs.wo = Transform4f(dr::rotate<Matrix4f>(Vector3f(0.f, 1.f, 0.f), theta_i))*Transform4f(dr::rotate<Matrix4f>(si.n, phi_i))*si.n;
        // bs.wo = si.wi;
        bs.pdf = 1.f;
        bs.eta = 1.f;
        bs.sampled_type = +BSDFFlags::DeltaReflection;
        bs.sampled_component = 0;

        active &= dr::neq(bs.pdf, 0.f) && Frame3f::cos_theta(bs.wo) > 0.f;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);

        return { bs, depolarizer<Spectrum>(value) & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        return 0.f;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        return 0.f;
    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DeltaReflection))
            return { 0.f, 0.f };

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        UnpolarizedSpectrum value =
            m_reflectance->eval(si, active);

        Float pdf = 1.f;

        return { depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f) };
    }

    Spectrum eval_diffuse_reflectance(const SurfaceInteraction3f &si,
                                      Mask active) const override {
        return m_reflectance->eval(si, active);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "naiveRetroReflector[" << std::endl
            << "  reflectance = " << string::indent(m_reflectance) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
    ref<Texture> m_alpha; // 拡がり角の正規分布の標準偏差
};

MI_IMPLEMENT_CLASS_VARIANT(naiveRetroReflector, BSDF)
MI_EXPORT_PLUGIN(naiveRetroReflector, "Naive RetroReflector material")
NAMESPACE_END(mitsuba)