#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/sampler.h>

NAMESPACE_BEGIN(mitsuba)
// diffuse BSDFをベースにBSSRDF化できないか試してみる。
template <typename Float, typename Spectrum>
class BSSRDFtest : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    BSSRDFtest(const Properties &props) : Base(props) {
        m_reflectance = props.texture<Texture>("reflectance", .5f);
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        deltaX = props.texture<Texture>("dx", 0.f)->mean();
        deltaY = props.texture<Texture>("dy", 0.f)->mean();
        deltaZ = props.texture<Texture>("dz", 0.f)->mean();
        randX = props.texture<Texture>("rx", 0.f)->mean();
        randY = props.texture<Texture>("ry", 0.f)->mean();
        randZ = props.texture<Texture>("rz", 0.f)->mean();
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("reflectance", m_reflectance.get(), +ParamFlags::Differentiable);
    }

    mitsuba::PCG32<UInt32> setRandomGenerator(const Float seed) const {
        mitsuba::PCG32<UInt32> rng(PCG32_DEFAULT_STATE, seed);
        rng.state = seed*1000000;
        rng.template next_float<Float>();
        return rng;
    }
    Float rand(mitsuba::PCG32<UInt32> &generator) const {
        return generator.template next_float<Float>();
    }
    Point3f rand3(mitsuba::PCG32<UInt32> &generator) const {
        Point3f sample3;
        sample3.x() = rand(generator);
        sample3.y() = rand(generator);
        sample3.z() = rand(generator);
        return sample3;
    }

    Float mag(Normal3f n) const {
        return dr::sqrt(n.x()*n.x() + n.y()*n.y() + n.z()*n.z());
    }
    Vector3f rotateVector(const Vector3f v, const Normal3f from, const Normal3f to) const {
        Normal3f fromNorm = normalize(from);
        Normal3f toNorm = normalize(to);
        Normal3f rotateAxis = normalize(cross(fromNorm, toNorm));
        Float s = mag(cross(fromNorm, toNorm));
        Float c = dot(fromNorm, toNorm);
        return c*v + (1-c)*dot(v,rotateAxis)*rotateAxis +s*cross(rotateAxis,v);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             SurfaceInteraction3f &si,
                                             Float sample1,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Point3f r3 = Point3f(0.f, 0.f, 0.f);
        auto rng = setRandomGenerator(sample1*1000000); // random number generator

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
            return { bs, 0.f };
        // printf("<%f, %f, %f>\t", si.p.x(), si.p.y(), si.p.z());
        // printf("<%f, %f, %f>\t", si.n.x(), si.n.y(), si.n.z());

        
        Vector3f c = Vector3f(deltaX, deltaY, deltaZ);
        r3 = rand3(rng);
        Vector3f randShift_local = Vector3f(randX*(r3.x()-0.5), randY*(r3.y()), randZ*(r3.z()-0.5));
        Vector3f randShift = rotateVector(randShift_local, Normal3f(0.f, 0.f, 1.f), si.n);
        si.p += c + randShift;
        bs.wo = reflect(si.wi);
        bs.pdf = 1.f;
        bs.eta = 1.f;
        bs.sampled_type = +BSDFFlags::DiffuseReflection;
        bs.sampled_component = 0;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);

        return { bs, depolarizer<Spectrum>(value) & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        return 0.f;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        return 0.f;
    }

    Spectrum eval_diffuse_reflectance(const SurfaceInteraction3f &si,
                                      Mask active) const override {
        return m_reflectance->eval(si, active);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "BSSRDFtest[" << std::endl
            << "  reflectance = " << string::indent(m_reflectance) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
    Float deltaX, deltaY, deltaZ;
    Float randX, randY, randZ;
};

MI_IMPLEMENT_CLASS_VARIANT(BSSRDFtest, BSDF)
MI_EXPORT_PLUGIN(BSSRDFtest, "BSSRDFtest")
NAMESPACE_END(mitsuba)
