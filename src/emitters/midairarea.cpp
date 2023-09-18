#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _emitter-midairarea:

Mid-air Area light (:monosp:`midairarea`)
---------------------------

.. pluginparameters::

 * - radiance
   - |spectrum| or |texture|
   - Specifies the emitted radiance in units of power per unit area per unit steradian.
   - |exposed|, |differentiable|

This plugin implements an mid-air area light, i.e. a light source that emits
diffuse illumination from the exterior of an arbitrary shape.
Since the emission profile of an area light is completely diffuse, it
has the same apparent brightness regardless of the observer's viewing
direction. Furthermore, since it occupies a nonzero amount of space, an
area light generally causes scene objects to cast soft shadows.

[NOTE] The luminance will attenuate according to the direction of incident ray 
       at the same rate as the luminance attenuation of mid-air image.

To create an area light source, simply instantiate the desired
emitter shape and specify an :monosp:`area` instance as its child:

.. tabs::
    .. code-tab:: xml
        :name: sphere-light

        <shape type="sphere">
            <emitter type="midairarea">
                <rgb name="radiance" value="1.0"/>
            </emitter>
        </shape>

    .. code-tab:: python

        'type': 'sphere',
        'emitter': {
            'type': 'area',
            'radiance': {
                'type': 'rgb',
                'value': 1.0,
            }
        }

 */

// Attenuation rates in angles from 10 to 60 degrees with 5 degree intervals
static constexpr float att_tbl[] = 
    {
        0.f, 0.07058824f, 0.12156863f, 0.17254902f, 0.24705882f, 0.34117647f, 
        0.42745098f, 0.52941176f, 0.56078431f, 0.44705882f, 0.33464567f
    };

template <typename Float, typename Spectrum>
class MidAirAreaLight final : public Emitter<Float, Spectrum> {
public:
    MI_IMPORT_BASE(Emitter, m_flags, m_shape, m_medium)
    MI_IMPORT_TYPES(Scene, Shape, Texture)

    using FloatStorage = DynamicBuffer<Float>;

    MidAirAreaLight(const Properties &props) : Base(props) {
        if (props.has_property("to_world"))
            Throw("Found a 'to_world' transformation -- this is not allowed. "
                  "The area light inherits this transformation from its parent "
                  "shape.");

        m_radiance = props.texture_d65<Texture>("radiance", 1.f);

        m_attenuation_table = dr::load<FloatStorage>(att_tbl, 10);

        m_flags = +EmitterFlags::Surface;
        if (m_radiance->is_spatially_varying())
            m_flags |= +EmitterFlags::SpatiallyVarying;
        dr::set_attr(this, "flags", m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        Base::traverse(callback);
        callback->put_object("radiance", m_radiance.get(), +ParamFlags::Differentiable);
    }

    Spectrum eval(const SurfaceInteraction3f &si, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);

        // Calculate the angle between the incident direction and the normal
        Float dp = Frame3f::cos_theta(si.wi);
        Float degree = dr::acos(dp) * 180 / dr::Pi<Float>;

        // Ignore the attenuation when the angle is out of range
        degree += (ScalarFloat) 45.0f;

        UInt32 i0 = dr::clamp(dr::floor(degree / 5), UInt32(0), UInt32(8));
        UInt32 i1 = i0 + 1;

        Float deg0 = Float (i0 * 5) + 10.0f;
        Float deg1 = Float (i1 * 5) + 10.0f;

        // Select attenuation rate according to the angle
        Float att0 = (Float)dr::gather<Float32>(m_attenuation_table, i0);
        Float att1 = (Float)dr::gather<Float32>(m_attenuation_table, i1);
        Float t           = (degree - deg0) / (deg1 - deg0);
        Float attenuation = (1 - t) * att0 + t * att1; // Linear interpolation
        //attenuation *= is_multi_coeff;
        
        // Calculate the radiance with attenuation according to the angle
        auto result = depolarizer<Spectrum>(m_radiance->eval(si, active) * attenuation) &
                      (Frame3f::cos_theta(si.wi) > 0.f);

        return result;
    }

    std::pair<Ray3f, Spectrum> sample_ray(Float time, Float wavelength_sample,
                                          const Point2f &sample2, const Point2f &sample3,
                                          Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleRay, active);

        // 1. Sample spatial component
        auto [ps, pos_weight] = sample_position(time, sample2, active);

        // 2. Sample directional component
        Vector3f local = warp::square_to_cosine_hemisphere(sample3);

        // 3. Sample spectral component
        SurfaceInteraction3f si(ps, dr::zeros<Wavelength>());
        auto [wavelength, wav_weight] =
            sample_wavelengths(si, wavelength_sample, active);
        si.time = time;
        si.wavelengths = wavelength;

        // Note: some terms cancelled out with `warp::square_to_cosine_hemisphere_pdf`.
        Spectrum weight = pos_weight * wav_weight * dr::Pi<ScalarFloat>;

        return { si.spawn_ray(si.to_world(local)),
                 depolarizer<Spectrum>(weight) };
    }

    std::pair<DirectionSample3f, Spectrum>
    sample_direction(const Interaction3f &it, const Point2f &sample, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSampleDirection, active);
        Assert(m_shape, "Can't sample from an area emitter without an associated Shape.");
        DirectionSample3f ds;
        SurfaceInteraction3f si;

        // One of two very different strategies is used depending on 'm_radiance'
        if (likely(!m_radiance->is_spatially_varying())) {
            // Texture is uniform, try to importance sample the shape wrt. solid angle at 'it'
            ds = m_shape->sample_direction(it, sample, active);
            active &= dr::dot(ds.d, ds.n) < 0.f && dr::neq(ds.pdf, 0.f);

            si = SurfaceInteraction3f(ds, it.wavelengths);
        } else {
            // Importance sample the texture, then map onto the shape
            auto [uv, pdf] = m_radiance->sample_position(sample, active);
            active &= dr::neq(pdf, 0.f);

            si = m_shape->eval_parameterization(uv, +RayFlags::All, active);
            si.wavelengths = it.wavelengths;
            active &= si.is_valid();

            ds.p = si.p;
            ds.n = si.n;
            ds.uv = si.uv;
            ds.time = it.time;
            ds.delta = false;
            ds.d = ds.p - it.p;

            Float dist_squared = dr::squared_norm(ds.d);
            ds.dist = dr::sqrt(dist_squared);
            ds.d /= ds.dist;

            Float dp = dr::dot(ds.d, ds.n);
            active &= dp < 0.f;
            ds.pdf = dr::select(active, pdf / dr::norm(dr::cross(si.dp_du, si.dp_dv)) *
                                        dist_squared / -dp, 0.f);
        }

        UnpolarizedSpectrum spec = m_radiance->eval(si, active) / ds.pdf;
        ds.emitter = this;
        return { ds, depolarizer<Spectrum>(spec) & active };
    }

    Float pdf_direction(const Interaction3f &it, const DirectionSample3f &ds,
                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);
        Float dp = dr::dot(ds.d, ds.n);
        active &= dp < 0.f;

        Float value;
        if (!m_radiance->is_spatially_varying()) {
            value = m_shape->pdf_direction(it, ds, active);
        } else {
            // This surface intersection would be nice to avoid..
            SurfaceInteraction3f si = m_shape->eval_parameterization(ds.uv, +RayFlags::dPdUV, active);
            active &= si.is_valid();

            value = m_radiance->pdf_position(ds.uv, active) * dr::sqr(ds.dist) /
                    (dr::norm(dr::cross(si.dp_du, si.dp_dv)) * -dp);
        }

        return dr::select(active, value, 0.f);
    }

    Spectrum eval_direction(const Interaction3f &it,
                            const DirectionSample3f &ds,
                            Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointEvaluate, active);
        Float dp = dr::dot(ds.d, ds.n);
        active &= dp < 0.f;

        SurfaceInteraction3f si(ds, it.wavelengths);
        UnpolarizedSpectrum spec = m_radiance->eval(si, active);
        return dr::select(active, depolarizer<Spectrum>(spec), 0.f);
    }

    std::pair<PositionSample3f, Float>
    sample_position(Float time, const Point2f &sample,
                    Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::EndpointSamplePosition, active);
        Assert(m_shape, "Cannot sample from an area emitter without an associated Shape.");

        // Two strategies to sample the spatial component based on 'm_radiance'
        PositionSample3f ps;
        if (!m_radiance->is_spatially_varying()) {
            // Radiance not spatially varying, use area-based sampling of shape
            ps = m_shape->sample_position(time, sample, active);
        } else {
            // Importance sample texture
            auto [uv, pdf] = m_radiance->sample_position(sample, active);
            active &= dr::neq(pdf, 0.f);

            auto si = m_shape->eval_parameterization(uv, +RayFlags::All, active);
            active &= si.is_valid();
            pdf /= dr::norm(dr::cross(si.dp_du, si.dp_dv));

            ps = si;
            ps.pdf = pdf;
            ps.delta = false;
        }

        Float weight = dr::select(active && ps.pdf > 0.f, dr::rcp(ps.pdf), Float(0.f));
        return { ps, weight };
    }

    std::pair<Wavelength, Spectrum>
    sample_wavelengths(const SurfaceInteraction3f &si, Float sample,
                       Mask active) const override {
        return m_radiance->sample_spectrum(
            si, math::sample_shifted<Wavelength>(sample), active);
    }

    ScalarBoundingBox3f bbox() const override { return m_shape->bbox(); }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "MidAirAreaLight[" << std::endl
            << "  radiance = " << string::indent(m_radiance) << "," << std::endl
            << "  surface_area = ";
        if (m_shape) oss << m_shape->surface_area();
        else         oss << "  <no shape attached!>";
        oss << "," << std::endl;
        if (m_medium) oss << string::indent(m_medium);
        else         oss << "  <no medium attached!>";
        oss << std::endl << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    ref<Texture> m_radiance;
    FloatStorage m_attenuation_table;
};

MI_IMPLEMENT_CLASS_VARIANT(MidAirAreaLight, Emitter)
MI_EXPORT_PLUGIN(MidAirAreaLight, "Mid-air Area emitter")
NAMESPACE_END(mitsuba)
