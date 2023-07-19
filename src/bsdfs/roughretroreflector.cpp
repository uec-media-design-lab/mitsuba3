#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/texture.h>


NAMESPACE_BEGIN(mitsuba)

/**!

.. _bsdf-null:

Null material (:monosp:`null`)
-------------------------------------------

This plugin models a completely invisible surface material.
Light will not interact with this BSDF in any way.

Internally, this is implemented as a forward-facing Dirac delta distribution.
Note that the standard :ref:`path tracer <integrator-path>` does not have a good sampling strategy to deal with this,
but the (:ref:`volumetric path tracer <integrator-volpath>`) does.

The main purpose of this material is to be used as the BSDF of a shape enclosing a participating medium.

 */

template <typename Float, typename Spectrum>
class RoughRetroreflector final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    RoughRetroreflector(const Properties &props) : Base(props) {
        m_components.push_back(BSDFFlags::Null | BSDFFlags::FrontSide | BSDFFlags::BackSide);
        m_flags = m_components.back();
        dr::set_attr(this, "flags", m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        if (!has_flag(m_flags, BSDFFlags::Anisotropic))
            callback->put_object("alpha", m_alpha_u_surface.get());
        else {
            callback->put_object("alpha_u_surface", m_alpha_u_surface.get());
            callback->put_object("alpha_v_surface", m_alpha_v_surface.get());
            callback->put_object("alpha_u_internal", m_alpha_u_internal.get());
            callback->put_object("alpha_v_internal", m_alpha_v_internal.get());
        }
        callback->put_parameter("eta", m_eta);
        if (m_surface_reflectance)
            callback->put_object("surface_reflectance", m_surface_reflectance.get());
        if (m_internal_reflectance)
            callback->put_object("internal_reflectance", m_internal_reflectance.get());
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /*sample1*/,
                                             const Point2f & /*sample2*/,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
        bool sample_transmission = ctx.is_enabled(BSDFFlags::Null, 0);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Spectrum result(0.f);
        if (sample_transmission) {
            bs.wo                = -si.wi;
            bs.sampled_component = 0;
            bs.sampled_type      = UInt32(+BSDFFlags::Null);
            bs.eta               = 1.f;
            bs.pdf               = 1.f;

            /* In an ordinary BSDF we would use depolarizer<Spectrum>(1.f) here
               to construct a depolarizing Mueller matrix. However, the null
               BSDF should leave the polarization state unaffected, and hence
               this is one of the few places where it is safe to directly use a
               scalar (which will broadcast to the identity matrix in polarized
               rendering modes). */
            result               = 1.f;
        }

        return { bs, result };
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

    Spectrum eval_null_transmission(const SurfaceInteraction3f & /*si*/,
                                    Mask /*active*/) const override {
        /* As above, we do not want the polarization state to change. So it is
           safe to return a scalar (which will broadcast to the identity
           matrix). */
        return 1.f;
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "RoughRetroeflector[" << std::endl
            << "  distribution = "     << m_type << "," << std::endl;

        oss << "  m_pa = " << m_pa << "," << std::endl
            << "  m_qa = " << m_qa << "," << std::endl
            << "  m_ra = " << m_ra << "," << std::endl
            << "  m_pb = " << m_pb << "," << std::endl
            << "  m_qb = " << m_qb << "," << std::endl
            << "  m_rb = " << m_rb << "," << std::endl;

        if (!has_flag(m_flags, BSDFFlags::Anisotropic)) {
            oss << "  alpha = "             << string::indent(m_alpha_v_surface) << "," << std::endl;
        } else {
            oss << "  alpha_u_surface = "   << string::indent(m_alpha_u_surface) << "," << std::endl
                << "  alpha_v_surface = "   << string::indent(m_alpha_v_surface) << "," << std::endl
                << "  alpha_u_internal = "  << string::indent(m_alpha_u_internal) << "," << std::endl
                << "  alpha_v_internal = "  << string::indent(m_alpha_v_internal) << "," << std::endl;
        }

        oss << "  surface_reflectance = "   << string::indent(m_surface_reflectance) << "," << std::endl;

        oss << "  internal_reflectance = " << string::indent(m_internal_reflectance) << ", " << std::endl;

        oss << "  eta = "       << m_eta << std::endl
            << "  eta_air = "   << m_eta_air << std::endl
            << "  eta_mat = "   << m_eta_mat << std::endl
            << "]";
        return oss.str();
    }
    MI_DECLARE_CLASS()
private:
    Normal3f m_pa, m_qa, m_ra;
    Normal3f m_pb, m_qb, m_rb;
    Normal3f n;
    ScalarFloat m_eta, m_inv_eta;
    ScalarFloat m_eta_mat, m_eta_air;
    ref<Texture> m_eta_base, m_k_base;
    ref<Texture> m_alpha_u_surface,  m_alpha_v_surface;
    ref<Texture> m_alpha_u_internal, m_alpha_v_internal;
    MicrofacetType m_type;
    bool m_sample_visible;
    Float diffuseFactor;
    Float m_surface_reflectance;
    Float m_internal_reflectance;
};

MI_IMPLEMENT_CLASS_VARIANT(RoughRetroreflector, BSDF)
MI_EXPORT_PLUGIN(RoughRetroreflector, "Rough retroreflector")
NAMESPACE_END(mitsuba)
