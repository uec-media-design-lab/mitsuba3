#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/texture.h>


NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class Halfmirror final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)
    Halfmirror(const Properties &props) : Base(props) {
        m_reflectance_texture   = props.texture<Texture>("reflectance_texture", 1.f);
        m_transmittance_texture   = props.texture<Texture>("transmittance_texture", 1.f);
        m_reflectance     = props.texture<Texture>("reflectance", 0.5f)->mean();
        m_transmittance   = props.texture<Texture>("transmittance", 0.5f)->mean();

        m_components.push_back(BSDFFlags::DeltaReflection | BSDFFlags::FrontSide |
                               BSDFFlags::BackSide);
        m_components.push_back(BSDFFlags::DeltaTransmission | BSDFFlags::FrontSide |
                               BSDFFlags::BackSide | BSDFFlags::NonSymmetric);

        m_flags = m_components[0] | m_components[1];
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_parameter("reflectance",   m_reflectance,   +ParamFlags::Differentiable);
        callback->put_parameter("transmittance", m_transmittance, +ParamFlags::Differentiable);

        if (m_reflectance_texture)
            callback->put_object("reflectance_texture",   m_reflectance_texture.get(),   +ParamFlags::Differentiable);
        if (m_transmittance_texture)
            callback->put_object("transmittance_texture", m_transmittance_texture.get(), +ParamFlags::Differentiable);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             SurfaceInteraction3f &si,
                                             Float sample1,
                                             const Point2f & /*sample2*/,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
        bool has_reflection   = ctx.is_enabled(BSDFFlags::DeltaReflection, 0),
             has_transmission = ctx.is_enabled(BSDFFlags::DeltaTransmission, 1);
        
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        UnpolarizedSpectrum weight = 0.f;

        Mask selected_r;
        Mask selected_t;

        if (likely(has_reflection && has_transmission)) {
            selected_r = (sample1 <= m_reflectance) && active;
            selected_t = !selected_r && (sample1 <= (m_reflectance + m_transmittance)) && active;

            bs.pdf = dr::select(selected_r, m_reflectance, m_transmittance);
        } else if (has_reflection || has_transmission) {
            selected_r = Mask(has_reflection) && active;
            selected_t = !selected_r && active;
            bs.pdf = 1.f;
        } else {
            return {bs, 0.f};
        }
        bs.sampled_component = dr::select(selected_r, UInt32(0), UInt32(1));
        bs.sampled_type      = dr::select(selected_r, UInt32(+BSDFFlags::DeltaReflection),
                                                      UInt32(+BSDFFlags::DeltaTransmission));
        // bs
        bs.wo = dr::select(selected_r, reflect(si.wi), bs.wo);
        bs.wo = dr::select(selected_t, -si.wi, bs.wo);
        bs.eta = 1.f;
        bs.sampled_component = dr::select(selected_r, UInt32(0), bs.sampled_component);
        bs.sampled_component = dr::select(selected_t, UInt32(1), bs.sampled_component);

        weight = dr::select(selected_r, m_reflectance * m_reflectance_texture->eval(si, selected_r), weight);
        weight = dr::select(selected_t, m_transmittance * m_transmittance_texture->eval(si, selected_t), weight);
        weight = 1.f;
        
        return { bs, weight };
    }

    Spectrum eval(const BSDFContext & /*ctx*/,
                  const SurfaceInteraction3f & /*si*/, const Vector3f & /*wo*/,
                  Mask /*active*/) const override {
        return 0.f;
    }

    Float pdf(const BSDFContext & /*ctx*/, const SurfaceInteraction3f & /*si*/,
              const Vector3f & /*wo*/, Mask /*active*/) const override {
        return 0.f;
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Half-mirror[" << std::endl;
        oss << "  reflectance = " << string::indent(m_reflectance) << "," << std::endl;
        oss << "  transmittance = " << string::indent(m_transmittance) << ", " << std::endl;
        if (m_reflectance_texture)
            oss << "  reflectance_texture = " << string::indent(m_reflectance_texture) << "," << std::endl;
        if (m_transmittance_texture)
            oss << "  transmittance_texture = " << string::indent(m_transmittance_texture) << ", " << std::endl;
        oss << "]";
        return oss.str();
    }
    MI_DECLARE_CLASS()
private:
    Float m_reflectance;
    Float m_transmittance;
    ref<Texture> m_reflectance_texture;
    ref<Texture> m_transmittance_texture;
};

MI_IMPLEMENT_CLASS_VARIANT(Halfmirror, BSDF)
MI_EXPORT_PLUGIN(Halfmirror, "Half-mirror material")
NAMESPACE_END(mitsuba)
