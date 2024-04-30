#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class LearnedMMAP final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    LearnedMMAP(const Properties &props) : Base(props) {
        m_flags = BSDFFlags::DeltaReflection | BSDFFlags::FrontSide |
                  BSDFFlags::BackSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);

        std::string material = props.string("material", "none");
        if (props.has_property("eta") || material == "none") {
            m_eta = props.texture<Texture>("eta", 1.0f);
            if (material != "none")
                Throw("Should specify eta or material not both.");
        } else {
            m_eta = props.texture<Texture>("eta", 1.0f);
        }
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("eta", m_eta.get(), +ParamFlags::Differentiable);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* sample1 */,
                                             const Point2f & /* sample2 */,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta > 0.0f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Spectrum value(0.0f);

        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::DeltaTransmission)))
            return { bs, value };

        bs.sampled_component = 0;
        bs.sampled_type      = +BSDFFlags::DeltaTransmission;
        // TODO: To be determined by learned smapling
        bs.wo  = retro_transmit(si.wi);
        bs.eta = 1.0f;
        bs.pdf = 1.0f;

        UnpolarizedSpectrum eta = m_eta->eval(si, active);
        UnpolarizedSpectrum retro_transmittance =
            m_retro_transmittance->eval(si, active);

        if constexpr (is_polarized_v<Spectrum>) {
            Vector3f wi_hat =
                         ctx.mode == TransportMode::Radiance ? bs.wo : si.wi,
                     wo_hat =
                         ctx.mode == TransportMode::Radiance ? si.wi : bs.wo;

            // Mueller matrix for specular reflection
            value = mueller::specular_reflection(
                UnpolarizedSpectrum(Frame3f::cos_theta(wi_hat)), eta);
            value = mueller::reverse(value);

            Vector3f n(0, 0, 1);
            Vector3f s_axis_in  = normalize(dr::cross(n, -wo_hat));
            Vector3f s_axis_out = normalize(dr::cross(n, wi_hat));

            value = mueller::rotate_mueller_basis(
                value, -wo_hat, s_axis_in, mueller::stokes_basis(-wo_hat),
                wi_hat, s_axis_out, mueller::stokes_basis(wi_hat));
            value *= mueller::absorber(retro_transmittance);
        } else {
            value = retro_transmittance;
        }

        return { bs, value & active };
    }

    Spectrum eval(const BSDFContext & /* ctx */,
                  const SurfaceInteraction3f & /* si */,
                  const Vector3f & /* wo */, Mask /* active */) const override {
        return 0.0f;
    }

    Float pdf(const BSDFContext & /* ctx */,
              const SurfaceInteraction3f & /* si */, const Vector3f & /* wo */,
              Mask /* active */) const override {
        return 0.0f;
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "LearnedMMAP[" << std::endl
            << "  eta = " << string::indent(m_reflectance) << std::endl
            << "  retro_transmittance = "
            << string::indent(m_retro_transmittance) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

private:
    ref<Texture> m_retro_transmittance;
    ref<Texture> m_eta;
};

MI_IMPLEMENT_CLASS_VARIANT(LearnedMMAP, BSDF)
MI_EXPORT_PLUGIN(LearnedMMAP, "Learned micro-mirror array plate material")
NAMESPACE_END(mitsuba)