#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/microfacet.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/traits.h>


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

 /* constants */
template <typename T> constexpr auto SqrtThree       = dr::scalar_t<T>(1.73205080756887729353);
template <typename T> constexpr auto InvSqrtThree    = dr::scalar_t<T>(0.57735026918962576451);

template <typename T> constexpr auto SqrtSix         = dr::scalar_t<T>(2.44948974278317809820);
template <typename T> constexpr auto InvSqrtSix      = dr::scalar_t<T>(0.40824829046386301637);

template <typename Float, typename Spectrum>
class RoughRetroreflector final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)

    RoughRetroreflector(const Properties &props) : Base(props) {
        m_pa = Normal3f(2*InvSqrtSix<Float>, 0, InvSqrtThree<Float>);                         // →
        m_qa = Normal3f(-InvSqrtSix<Float>, dr::InvSqrtTwo<Float>, InvSqrtThree<Float>);    // ↖
        m_ra = Normal3f(-InvSqrtSix<Float>, -dr::InvSqrtTwo<Float>, InvSqrtThree<Float>);   // ↙
        m_pb = Normal3f(-2*InvSqrtSix<Float>, 0, InvSqrtThree<Float>);                       // ←
        m_qb = Normal3f(InvSqrtSix<Float>, dr::InvSqrtTwo<Float>, InvSqrtThree<Float>);    // ↗
        m_rb = Normal3f(InvSqrtSix<Float>, -dr::InvSqrtTwo<Float>, InvSqrtThree<Float>);   // ↘

        // # iors
        ScalarFloat int_ior = lookup_ior(props, "int_ior", "bk7");
        ScalarFloat ext_ior = lookup_ior(props, "ext_ior", "air");
        std::string material = props.string("base_material", "none");
        if (props.has_property("base_eta") || material == "none") {
            // etaで指定されてる場合
            // デフォでは屈折率実部0, 虚部1で完全に吸収する素材が指定される
            m_eta_base = props.texture<Texture>("base_eta", 0.f);
            m_k_base   = props.texture<Texture>("base_k",   1.f);
            if (material != "none")
                Throw("Should specify either (eta, k) or material, not both.");
        } else {
            // 材質が指定されてる場合は参照する。見つからなければCu
            std::tie(m_eta_base, m_k_base) = complex_ior_from_file<Spectrum, Texture>(props.string("base_material", "Al"));
        }

        if (int_ior < 0.f || ext_ior < 0.f || int_ior == ext_ior)
            Throw("The interior and exterior indices of "
                  "refraction must be positive and differ!");
        m_eta = int_ior / ext_ior;
        m_inv_eta = ext_ior / int_ior;
        m_eta_mat = int_ior;
        m_eta_air = ext_ior;

        // # distributions
        if (props.has_property("distribution")) {
            std::string distr = string::to_lower(props.string("distribution"));
            if (distr == "beckmann")
                m_type = MicrofacetType::Beckmann;
            else if (distr == "ggx")
                m_type = MicrofacetType::GGX;
            else
                Throw("Specified an invalid distribution \"%s\", must be "
                      "\"beckmann\" or \"ggx\"!", distr.c_str());
        } else {
            m_type = MicrofacetType::Beckmann;
        }

        if (props.has_property("alpha_u_surface") || props.has_property("alpha_v_surface")) {
            if (!props.has_property("alpha_u_surface") || !props.has_property("alpha_v_surface"))
                Throw("Microfacet model: both 'alpha_u' and 'alpha_v' must be specified.");
            if (props.has_property("alpha_surface"))
                Throw("Microfacet model: please specify"
                      "either 'alpha' or 'alpha_u_surface'/'alpha_v_surface'.");
            m_alpha_u_surface = props.texture<Texture>("alpha_u_surface");
            m_alpha_v_surface = props.texture<Texture>("alpha_v_surface");
        } else if (props.has_property("alpha_surface")) {
            m_alpha_u_surface = m_alpha_v_surface = props.texture<Texture>("alpha_surface", 0.1f);
        }
        if (props.has_property("alpha_u_internal") || props.has_property("alpha_v_internal")) {
            if (!props.has_property("alpha_u_internal") || !props.has_property("alpha_v_internal"))
                Throw("Microfacet model: both 'alpha_u' and 'alpha_v' must be specified.");
            if (props.has_property("alpha_internal"))
                Throw("Microfacet model: please specify"
                      "either 'alpha' or 'alpha_u_internal'/'alpha_v_internal'.");
            m_alpha_u_internal = props.texture<Texture>("alpha_u_internal");
            m_alpha_v_internal = props.texture<Texture>("alpha_v_internal");
        } else if (props.has_property("alpha_internal")) {
            m_alpha_u_internal = m_alpha_v_internal = props.texture<Texture>("alpha_internal", 0.1f);
        }
        if (props.has_property("alpha")) {
            m_alpha_u_surface = m_alpha_v_surface = m_alpha_u_internal = m_alpha_v_internal =props.texture<Texture>("alpha", 0.1f);
        }

        // # flags
        BSDFFlags anistropic = ((m_alpha_u_surface != m_alpha_v_surface) || (m_alpha_u_internal != m_alpha_v_internal)) ? BSDFFlags::Anisotropic : BSDFFlags(0);
        // # reflection, retroreflection, diffuse
        m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide | BSDFFlags::BackSide | anistropic);
        m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::GlossyTransmission | BSDFFlags::FrontSide | BSDFFlags::BackSide | anistropic);
        m_components.push_back(BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide);
        m_flags= m_components[0] | m_components[1] | m_components[2];
        dr::set_attr(this, "flags", m_flags);
    
        // sample_visible
        m_sample_visible = props.get<bool>("sample_visible", true);

        // diffuseFactor
        Float F_d = 0.f;
        if (m_eta<1.f) {
            F_d = 0.919317 - 3.4793*m_eta + 6.75335*m_eta*m_eta - 7.80989*m_eta*m_eta*m_eta + 4.98554*m_eta*m_eta*m_eta*m_eta - 1.36881*m_eta*m_eta*m_eta*m_eta*m_eta;
        } else {
            F_d = -9.23372 + 22.2272*m_eta - 20.9292*m_eta*m_eta + 10.2291*m_eta*m_eta*m_eta - 2.54396*m_eta*m_eta*m_eta*m_eta + 0.254913*m_eta*m_eta*m_eta*m_eta*m_eta;
        }
        diffuseFactor = dr::InvPi<Float> * dr::sqr(m_inv_eta) / (1.f - F_d);
        
        // surface normal(ideal)
        n = Normal3f(0.f, 0.f, 1.f);

        // reflectance
        m_surface_reflectance   = props.texture<Texture>("surface_reflectance", 1.f);
    
        m_internal_reflectance   = props.texture<Texture>("internal_reflectance", 1.f);
        randX = props.texture<Texture>("rx", 0.f)->mean();
        randY = props.texture<Texture>("ry", 0.f)->mean();
        randZ = props.texture<Texture>("rz", 0.f)->mean();
    }

    void traverse(TraversalCallback *callback) override {
        if (!has_flag(m_flags, BSDFFlags::Anisotropic))
            callback->put_object("alpha", m_alpha_u_surface.get(),                  ParamFlags::Differentiable | ParamFlags::Discontinuous);
        else {
            callback->put_object("alpha_u_surface", m_alpha_u_surface.get(),        ParamFlags::Differentiable | ParamFlags::Discontinuous);
            callback->put_object("alpha_v_surface", m_alpha_v_surface.get(),        ParamFlags::Differentiable | ParamFlags::Discontinuous);
            callback->put_object("alpha_u_internal", m_alpha_u_internal.get(),      ParamFlags::Differentiable | ParamFlags::Discontinuous);
            callback->put_object("alpha_v_internal", m_alpha_v_internal.get(),      ParamFlags::Differentiable | ParamFlags::Discontinuous);
        }
        callback->put_parameter("eta", m_eta, ParamFlags::Differentiable | ParamFlags::Discontinuous);
        if (m_surface_reflectance)
            callback->put_object("surface_reflectance", m_surface_reflectance.get(),   +ParamFlags::Differentiable);
        if (m_internal_reflectance)
            callback->put_object("internal_reflectance", m_internal_reflectance.get(), +ParamFlags::Differentiable);
    }

    
    // print vector(or normal)(s)
    // void printV(const Vector3f v) const {
    //     printf("(%f, %f, %f)\t", v.x(), v.y(), v.z());
    // }
    // void printN(const Normal3f n) const {
    //     printf("<%f, %f, %f>\t", n.x(), n.y(), n.z());
    // }

    // random numbers
    mitsuba::PCG32<UInt32> setRandomGenerator(const Float seed) const {
        mitsuba::PCG32<UInt32> rng(PCG32_DEFAULT_STATE, seed);
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
    Point3f rand3c(mitsuba::PCG32<UInt32> &generator) const {
        Point3f sample3;
        sample3.x() = rand(generator)-0.5f;
        sample3.y() = rand(generator)-0.5f;
        sample3.z() = rand(generator)-0.5f;
        return sample3;
    }

    // rotations
    Float mag(Normal3f n) const {
        return sqrt(n.x()*n.x() + n.y()*n.y() + n.z()*n.z());
    }
    Normal3f rotateNormal(const Normal3f n, const Normal3f from, const Normal3f to) const {
        Normal3f fromNorm = normalize(from);
        Normal3f toNorm = normalize(to);
        Normal3f rotateAxis = normalize(cross(fromNorm, toNorm));
        Float s = mag(cross(fromNorm, toNorm));
        Float c = dot(fromNorm, toNorm);
        return c*n + (1-c)*dot(n,rotateAxis)*rotateAxis +s*cross(rotateAxis,n);
    }
    Vector3f rotateVector(const Vector3f v, const Normal3f from, const Normal3f to) const {
        Normal3f fromNorm = normalize(from);
        Normal3f toNorm = normalize(to);
        Normal3f rotateAxis = normalize(cross(fromNorm, toNorm));
        Float s = mag(cross(fromNorm, toNorm));
        Float c = dot(fromNorm, toNorm);
        return c*v + (1-c)*dot(v,rotateAxis)*rotateAxis +s*cross(rotateAxis,v);
    }

    // sampleRoute
    std::tuple<Normal3f, Normal3f, Normal3f, Float> sampleRoute(const Vector3f wi, const Normal3f a, const Normal3f b, const Normal3f c, mitsuba::PCG32<UInt32> &rng) const {
        // wi toward a,b,c normals
        // pa,qa,raもしくはpb,qb,rbが与えられたときにサンプルする順番とその確率を与える
        Float prob=0.5f;    // △ or ▽
        Normal3f n1, n2, n3;
        Normal3f nt1, nt2;
        Vector3f vt;
        Point2f r2 = rand2(rng);

        // 1st surface
        Float inner_a = dr::clamp(dot(-wi, a), 0.f, 1.f);
        Float inner_b = dr::clamp(dot(-wi, b), 0.f, 1.f);
        Float inner_c = dr::clamp(dot(-wi, c), 0.f, 1.f);
        Float invDenominator = 1.f / (inner_a + inner_b + inner_c);
        n1 =    dr::select(Mask(r2.x() <= inner_a*invDenominator), a,
                dr::select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), b, c));
        prob =  dr::select(Mask(r2.x() <= inner_a*invDenominator), inner_a,
                dr::select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), inner_b, inner_c)) * invDenominator * 0.5f;

        // 2nd, 3rd surfaces
        nt1 =   dr::select(Mask(r2.x() <= inner_a*invDenominator), b, 
                dr::select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), a, a));
        nt2 =   dr::select(Mask(r2.x() <= inner_a*invDenominator), c, 
                dr::select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), c, b));
        vt = reflect(-wi, n1);
        Float inner_nt1 = dr::clamp(dot(-vt, nt1), 0.f, 1.f);
        Float inner_nt2 = dr::clamp(dot(-vt, nt2), 0.f, 1.f);
        invDenominator = 1.f / (inner_nt1 + inner_nt2);

        n2 =            dr::select(Mask(r2.y() <= inner_nt1*invDenominator), nt1, nt2);
        n3 =            dr::select(Mask(r2.y() <= inner_nt1*invDenominator), nt2, nt1);
        prob = prob *   dr::select(Mask(r2.y() <= inner_nt1*invDenominator), inner_nt1, inner_nt2) * invDenominator;

        return std::forward_as_tuple(n1, n2, n3, prob);
    }

    // ERA
    Float ERA(const Float cos_theta_i, const Float sin_phi_i, const Float cos_phi_i) const {
        Float x = dr::rad_to_deg(acos(abs(cos_theta_i)));
        Float x_phi = dr::mulsign(dr::rad_to_deg(acos(cos_phi_i)), sin_phi_i) + 360.f;  // 0 -- 360 [deg]
        Float x_phi_blend = abs(fmod(x_phi, 60.f) - 30.f) * 0.1f * 0.33333333333333333333f;
        Float val_upper = ERA_upperBoundary(x);
        Float val_lower = ERA_lowerBoundary(x);
        // printf("phi:%f\tphi:%f\tphimod:%f\tphimodabs30:%f\tblend:%f\n",dr::rad_to_deg(acos(cos_phi_i)), x_phi, fmod(x_phi, 60.f), abs(fmod(x_phi, 60.f) - 30.f), x_phi_blend);
        // printf("cos:%f\tsin:%f\tphi:%f\tphi':%f\tphi-blend:%f\n", cos_phi_i, sin_phi_i, x_phi, dr::rad_to_deg(acos(cos_phi_i)), x_phi_blend);
        // printf("theta:%f\tphi':%f\tval:%f\n", x, x_phi_blend, x_phi_blend*val_upper + (1.f-x_phi_blend)*val_lower);
        return x_phi_blend*val_upper + (1.f-x_phi_blend)*val_lower;
    }

    Float ERA_lowerBoundary(const Float x) const {
        Float a = -0.480967;
        Float b = 0.085856;
        Float c = -2.860831;
        Float d = 0.084810;
        Float val = a*atan(b*x +c)+d;
        // printf("x:%f, val:%f\n", x, val);
        return dr::select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
    }
    Float ERA_upperBoundary(const Float x) const {
        Float a = -0.353394;
        Float b = 0.055244;
        Float c = -1.632606;
        Float d = 0.322462;
        Float val = a*atan(b*x +c)+d;
        // printf("x:%f, val:%f\n", x, val);
        return dr::select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
    }

    // class
    class Path {
        // 
        public:
            // MI_IMPORT_BASE(BSDF, m_flags, m_components)
            // MI_IMPORT_BASE(m_components)
            MI_IMPORT_TYPES(Texture, MicrofacetDistribution)
            // MI_IMPORT_TYPES(MicrofacetType)
            Path(Normal3f n, Normal3f n1, Normal3f n2, Normal3f n3, Vector3f wi, Vector3f wo, ScalarFloat eta_air, ScalarFloat eta_mat, ref<Texture> base_eta, ref<Texture> base_k) {
                this->n = n;
                this->ni = n;
                this->n1 = n1;
                this->n2 = n2;
                this->n3 = n3;
                this->no = n;
                this->ea = eta_air;
                this->em = eta_mat;
                this->reta = em / ea;
                this->invreta = ea / em;
                this->eta_base = base_eta;
                this->k_base = base_k;

                // 評価用の平均経路
                Vector3f wii, wv1i, wv2i, wv3i, wv4i, woi;
                Vector3f wio, wv1o, wv2o, wv3o, wv4o, woo;
                std::tie(wii, wv1i, wv2i, wv3i, wv4i, woi) = getPath(wi, ni, n1, n2, n3, no);
                std::tie(woo, wv4o, wv3o, wv2o, wv1o, wio) = getPath(wo, no, n3, n2, n1, ni);
                this->mui  = normalize(wii - wio);
                this->muv1 = normalize(wv1i - wv1o);
                this->muv2 = normalize(wv2i - wv2o);
                this->muv3 = normalize(wv3i - wv3o);
                this->muv4 = normalize(wv4i - wv4o);
                this->muo  = normalize(woi - woo);

                // 平均経路基準の法線ベクトル
                this->mi = normalize(mui*ea - muv1*em);
                this->m1 = normalize(-muv1 + muv2);
                this->m2 = normalize(-muv2 + muv3);
                this->m3 = normalize(-muv3 + muv4);
                this->mo = normalize(muv4*em - muo*ea);
                
                this->pathProb = getProb();
            }
            std::tuple<Vector3f, Vector3f, Vector3f, Vector3f, Vector3f, Vector3f> getPath(Vector3f i, Normal3f ni, Normal3f n1, Normal3f n2, Normal3f n3, Normal3f no) const {
                // 常に進行方向を向き続けるベクトルを返す(統制)
                Vector3f vi, v1, v2, v3, v4, vo;
                vi = -i;
                
                auto[r, cos_theta_t, eta_it, eta_ti] = fresnel(Frame3f::cos_theta(i), Float(reta));
                v1 = refract(i, cos_theta_t, eta_ti);
                v2 = reflect(-v1, n1);
                v3 = reflect(-v2, n2);
                v4 = reflect(-v3, n3);
                auto [ro, cos_theta_o, eta_ito, eta_tio] = fresnel(dot(-v4, n), Float(reta));
                vo = refract(-v4, n, cos_theta_o, eta_tio);
            
                return std::forward_as_tuple(vi, v1, v2, v3, v4, vo);
            }
            // class - rotations
            Float mag(Normal3f n) const {
                return sqrt(n.x()*n.x() + n.y()*n.y() + n.z()*n.z());
            }
            Normal3f rotateNormal(const Normal3f n, const Normal3f from, const Normal3f to) const {
                Normal3f fromNorm = normalize(from);
                Normal3f toNorm = normalize(to);
                Normal3f rotateAxis = normalize(cross(fromNorm, toNorm));
                Float s = mag(cross(fromNorm, toNorm));
                Float c = dot(fromNorm, toNorm);
                return c*n + (1-c)*dot(n,rotateAxis)*rotateAxis +s*cross(rotateAxis,n);
            }
            Vector3f rotateVector(const Vector3f v, const Normal3f from, const Normal3f to) const {
                Normal3f fromNorm = normalize(from);
                Normal3f toNorm = normalize(to);
                Normal3f rotateAxis = normalize(cross(fromNorm, toNorm));
                Float s = mag(cross(fromNorm, toNorm));
                Float c = dot(fromNorm, toNorm);
                return c*v + (1-c)*dot(v,rotateAxis)*rotateAxis +s*cross(rotateAxis,v);
            }

            // class - weight
            Float getProb() {
                Float prob = 0.5f;
                // n1,n2,n3を使っている。m1,m2,m3の方が良い？
                ////////////////////NAN check////////////////////////////////
                Mask valid = true;

                // 1st surface
                Float inner_1st_n1 = dr::clamp(dot(-muv1, m1), 0.f, 1.f);
                Float inner_1st_n2 = dr::clamp(dot(-muv1, m2), 0.f, 1.f);
                Float inner_1st_n3 = dr::clamp(dot(-muv1, m3), 0.f, 1.f);
                Float invDenominator1st = 1.f / (inner_1st_n1 + inner_1st_n2 + inner_1st_n3);

                prob = prob * inner_1st_n1*invDenominator1st;

                Float inner_2nd_n2 = dr::clamp(dot(-muv2, m2), 0.f, 1.f);
                Float inner_2nd_n3 = dr::clamp(dot(-muv2, m3), 0.f, 1.f);
                Float invDenominator2nd = 1.f / (inner_2nd_n2 + inner_2nd_n3);
                
                prob = prob * inner_2nd_n2*invDenominator2nd;
                
                valid = Mask(inner_1st_n1 + inner_1st_n2 + inner_1st_n3 > 0.f)&& valid;
                valid = Mask(inner_2nd_n2 + inner_2nd_n3>0.f)&& valid;
                valid = Mask(isfinite(prob)) && valid;
                return dr::select(valid, prob, 0.f);
            }

            // class - ERA
            Float ERA(const Float cos_theta_i, const Float sin_phi_i, const Float cos_phi_i) const {
            Float x = dr::rad_to_deg(acos(abs(cos_theta_i)));
            Float x_phi = dr::mulsign(dr::rad_to_deg(acos(cos_phi_i)), sin_phi_i) + 360.f;  // 0 -- 360 [deg]
            Float x_phi_blend = abs(fmod(x_phi, 60.f) - 30.f) * 0.1f * 0.33333333333333333333f;
            Float val_upper = ERA_upperBoundary(x);
            Float val_lower = ERA_lowerBoundary(x);
            return x_phi_blend*val_upper + (1.f-x_phi_blend)*val_lower;
        }

        Float ERA_lowerBoundary(const Float x) const {
            Float a = -0.480967;
            Float b = 0.085856;
            Float c = -2.860831;
            Float d = 0.084810;
            Float val = a*atan(b*x +c)+d;
            // printf("x:%f, val:%f\n", x, val);
            return dr::select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
        }
        Float ERA_upperBoundary(const Float x) const {
            Float a = -0.353394;
            Float b = 0.055244;
            Float c = -1.632606;
            Float d = 0.322462;
            Float val = a*atan(b*x +c)+d;
            // printf("x:%f, val:%f\n", x, val);
            return dr::select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
        }

        std::tuple<Spectrum, Spectrum> eval(const SurfaceInteraction3f &si, ref<Texture> alpha_surface_u, ref<Texture> alpha_surface_v, ref<Texture> alpha_internal_u, ref<Texture> alpha_internal_v, MicrofacetType mType, Mask active, bool sample_visible) const {
            Spectrum val_rr = 0.f;
            Spectrum val_d  = 0.f;

            MicrofacetDistribution distr_surface(mType, alpha_surface_u->eval_1(si, active), alpha_surface_v->eval_1(si, active), sample_visible);
            MicrofacetDistribution distr_internal(mType, alpha_internal_u->eval_1(si, active), alpha_internal_v->eval_1(si, active), sample_visible);

            // Walter's trick
            if (unlikely(!sample_visible)) {
                distr_surface.scale_alpha(1.2f - .2f * sqrt(abs(Frame3f::cos_theta(si.wi))));
                distr_internal.scale_alpha(1.2f - .2f * sqrt(abs(Frame3f::cos_theta(si.wi))));
            }

            // D
            Float Di = distr_surface.eval(mi);
            Float D1 = distr_internal.eval(rotateNormal(m1, n1, ni));
            Float D2 = distr_internal.eval(rotateNormal(m2, n2, ni));
            Float D3 = distr_internal.eval(rotateNormal(m3, n3, ni));
            Float Do = distr_surface.eval(mo);

            // F
            dr::Complex<UnpolarizedSpectrum> base_eta_k(eta_base->eval(si, active),
                                                    k_base->eval(si, active));
            Float Fi = std::get<0>(fresnel(dot(-mui, mi), Float(reta)));
            auto F1 = fresnel_conductor(UnpolarizedSpectrum(dot(-muv1, m1)), base_eta_k);
            auto F2 = fresnel_conductor(UnpolarizedSpectrum(dot(-muv2, m2)), base_eta_k);
            auto F3 = fresnel_conductor(UnpolarizedSpectrum(dot(-muv3, m3)), base_eta_k);
            Float Fo = std::get<0>(fresnel(dot(muo, mo), Float(reta)));
            // incorrect Path check
            // printf("(%f:  %f, %f, %f, %f, %f)\n", pathProb, Fi, F1, F2, F3, Fo);
            Mask incorrectRRpath = !isfinite(Fo);
            Fi = dr::select(incorrectRRpath, 0.f, Fi);
            F1 = dr::select(incorrectRRpath, 0.f, F1);
            F2 = dr::select(incorrectRRpath, 0.f, F2);
            F3 = dr::select(incorrectRRpath, 0.f, F3);
            Fo = dr::select(incorrectRRpath, 0.f, Fo);

            // G
            Float Gi = distr_surface.G(-mui, muv1, mi);
            Float G1 = distr_internal.G(rotateVector(-muv1,n1,ni), rotateVector(muv2,n1,ni), rotateNormal(m1,n1,ni));
            Float G2 = distr_internal.G(rotateVector(-muv2,n2,ni), rotateVector(muv3,n2,ni), rotateNormal(m2,n2,ni));
            Float G3 = distr_internal.G(rotateVector(-muv3,n3,ni), rotateVector(muv4,n3,ni), rotateNormal(m3,n3,ni));
            Float Go = distr_surface.G(-muv4, muo, mo);

            // J
            Float Ji = abs(ea*ea*dot(mui, mi))*abs(dot(muv1, mi)) / (abs(dot(mui, ni))*abs(dot(muv1, ni)) * dr::sqr(ea*dot(-mui,mi) + em*dot(muv1,mi)));
            Float J1 = 0.25f / (abs(dot(muv1,n1)*dot(muv2,n1)));
            Float J2 = 0.25f / (abs(dot(muv2,n2)*dot(muv3,n2)));
            Float J3 = 0.25f / (abs(dot(muv3,n3)*dot(muv4,n3)));
            Float Jo = abs(em*em*dot(muv4, mo))*abs(dot(muo, mo)) / (abs(dot(muv4, no))*abs(dot(muo, no)) * dr::sqr(em*dot(-muv4,mo) + ea*dot(muo,mo)));

            // ERA
            Float era_cos_i = Frame3f::cos_theta(-muv1);
            auto [era_sin_phi, era_cos_phi] = Frame3f::sincos_phi(-muv1);
            Float era = ERA(era_cos_i, era_sin_phi, era_cos_phi);

            // RetroReflection
            // val_rr = pathProb * (1.f-Fi)*F1*F2*F3*(1.f-Fo) * pow(Gi*G1*G2*G3*Go * Di*D1*D2*D3*Do * Ji*J1*J2*J3*Jo, 0.2);
            val_rr = dr::select(!incorrectRRpath, pathProb * (1.f-Fi)*F1*F2*F3*(1.f-Fo) * Gi*G1*G2*G3*Go * Di*D1*D2*D3*Do * Ji*J1*J2*J3*Jo * era, 0.f);
            // printf("[P:%f, F:%f, G:%f, D:%f, J:%f]\n", pathProb, (1.f-Fi)*F1*F2*F3*(1.f-Fo), Gi*G1*G2*G3*Go, pow(Di*D1*D2*D3*Do, 0.2), Ji*J1*J2*J3*Jo);

            // Diffuse
            // diffuseFactorを除いた部分
            val_d = dr::select(isfinite(pathProb), (1.f-Fi)* (era * pathProb * (1.f-F1*F2*F3*(1.f-Fo)) + (1.f-era)), 0.f);
            return std::forward_as_tuple(val_rr, val_d);
        }
        std::tuple<Float, Float> pdf(const SurfaceInteraction3f &si, ref<Texture> alpha_surface_u, ref<Texture> alpha_surface_v, ref<Texture> alpha_internal_u, ref<Texture> alpha_internal_v, MicrofacetType mType, Mask active, bool sample_visible) const {
            Float prob_rr = 0.f;
            Float prob_d  = 0.f;

            MicrofacetDistribution distr_surface(mType, alpha_surface_u->eval_1(si, active), alpha_surface_v->eval_1(si, active), sample_visible);
            MicrofacetDistribution distr_internal(mType, alpha_internal_u->eval_1(si, active), alpha_internal_v->eval_1(si, active), sample_visible);
            // Walter's trick
            if (unlikely(!sample_visible)) {
                distr_surface.scale_alpha(1.2f - .2f * sqrt(abs(Frame3f::cos_theta(si.wi))));
                distr_internal.scale_alpha(1.2f - .2f * sqrt(abs(Frame3f::cos_theta(si.wi))));
            }

            // D
            Float Di = distr_surface.pdf(abs(dot(mui, ni)), mi);
            Float D1 = distr_internal.pdf(abs(dot(muv1, n1)), rotateNormal(m1, n1, ni));
            Float D2 = distr_internal.pdf(abs(dot(muv2, n2)), rotateNormal(m2, n2, ni));
            Float D3 = distr_internal.pdf(abs(dot(muv3, n3)), rotateNormal(m3, n3, ni));
            Float Do = distr_surface.pdf(abs(dot(muv4, no)), mo);

            // F
            dr::Complex<UnpolarizedSpectrum> base_eta_k(eta_base->eval(si, active),
                                                    k_base->eval(si, active));
            Float Fi = std::get<0>(fresnel(dot(-mui, mi), Float(reta)));
            Float F1 = 1.f;
            Float F2 = 1.f;
            Float F3 = 1.f;
            Float Fo = std::get<0>(fresnel(dot(muo, mo), Float(reta)));
            // incorrect Path check
            Mask incorrectRRpath = !isfinite(Fo);
            F1 = dr::select(incorrectRRpath, 0.f, F1);
            F2 = dr::select(incorrectRRpath, 0.f, F2);
            F3 = dr::select(incorrectRRpath, 0.f, F3);
            Fo = dr::select(incorrectRRpath, 0.f, Fo);

            // J
            Float Ji = em*em*abs(dot(muv1, mi)) / dr::sqr(ea*dot(-mui,mi) + em*dot(muv1,mi));
            Float J1 = 0.25f / abs(dot(muv2,n1));
            Float J2 = 0.25f / abs(dot(muv3,n2));
            Float J3 = 0.25f / abs(dot(muv4,n3));
            Float Jo = ea*ea*abs(dot(muo, mo)) / dr::sqr(em*dot(-muv4,mo) + ea*dot(muo,mo));

            // ERA
            Float era_cos_i = Frame3f::cos_theta(-muv1);
            auto [era_sin_phi, era_cos_phi] = Frame3f::sincos_phi(-muv1);
            Float era = ERA(era_cos_i, era_sin_phi, era_cos_phi);

            // RetroReflection
            prob_rr = dr::select(!incorrectRRpath, 0.f, pathProb * (1.f-Fi)*F1*F2*F3*(1.f-Fo) * Di*D1*D2*D3*Do * Ji*J1*J2*J3*Jo * era);

            // Diffuse
            // diffuseFactorを除いた部分
            prob_d = (1.f-Fi)* (pathProb * (era * 1.f-F1*F2*F3*(1.f-Fo)) + (1.f-era)) * warp::square_to_cosine_hemisphere_pdf(muo);

            return std::forward_as_tuple(prob_rr, prob_d);
        }


        private:
        Normal3f n, ni, no, n1, n2, n3;
        Normal3f mi, mo, m1, m2, m3;
        Vector3f mui, muv1, muv2, muv3, muv4, muo;
        ScalarFloat ea, em;
        ScalarFloat reta, invreta;
        Float pathProb;
        ref<Texture> eta_base, k_base;
    };


    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             SurfaceInteraction3f &si,
                                             Float sample1,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);
        bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        bool has_retroreflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        has_retroreflection &= ctx.is_enabled(BSDFFlags::GlossyTransmission, 1);
        bool has_diffusereflection = ctx.is_enabled(BSDFFlags::DiffuseReflection, 2);
        // printf("asu: %f\tasv: %f\taiu: %f\taiv: %f\n",
        // m_alpha_u_surface->eval_1(si, active),
        // m_alpha_v_surface->eval_1(si, active),
        // m_alpha_u_internal->eval_1(si, active),
        // m_alpha_v_internal->eval_1(si, active));
        Float r1 = 0.f;                 // 汎用乱数(1D)
        Point2f r2 = Point2f(0.f, 0.f); // 汎用乱数(2D)

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        Spectrum weight = 0.f;

        auto rng = setRandomGenerator(sample1*1000000);
        MicrofacetDistribution distr_surface(m_type, m_alpha_u_surface->eval_1(si, active), m_alpha_v_surface->eval_1(si, active), m_sample_visible);
        MicrofacetDistribution sample_distr_surface(distr_surface);
        MicrofacetDistribution distr_internal(m_type, m_alpha_u_internal->eval_1(si, active), m_alpha_v_internal->eval_1(si, active), m_sample_visible);
        MicrofacetDistribution sample_distr_internal(distr_internal);
        
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        active &= cos_theta_i != 0.f;
        
        if (unlikely(!m_sample_visible)) {  // Walter's trick
            sample_distr_surface.scale_alpha(1.2f - .2f * sqrt(abs(cos_theta_i)));
            sample_distr_internal.scale_alpha(1.2f - .2f * sqrt(abs(cos_theta_i)));
        }

        // 最初の表面
        auto [mi, Di] = sample_distr_surface.sample(dr::mulsign(si.wi, cos_theta_i), sample2);
        auto [r_i, cos_theta_v1, eta_i_v1, eta_v1_i] = fresnel(dot(si.wi, mi), Float(m_eta));
        Float Fi = r_i;
        Vector3f wv1 = refract(si.wi, mi, cos_theta_v1, eta_v1_i);

        // 表面反射かそれ以外か分ける
        Mask selected_r = Mask(sample1 <= Fi) && active;
        Mask selected_t = Mask(sample1 >  Fi) && active;
        Mask selected_rr = false;   // decide later
        Mask selected_d = false;

        // 内部反射のルートを計算する
        Vector3f wo_rr = Vector3f(0.f, 0.f, -1.f);
        Float pathProb = 0.f;
        Spectrum F1 = 0.f;
        Spectrum F2 = 0.f;
        Spectrum F3 = 0.f;
        Float Fo = 0.f;
        Float D1 = 0.f;
        Float D2 = 0.f;
        Float D3 = 0.f;
        Float Do = 0.f;
        // Vectors, Normals
        Normal3f n1 = n, n2 = n, n3 = n;    // ideal
        Normal3f m1 = n, m2 = n, m3 = n, mo = n;    // sampled
        Vector3f wv2, wv3, wv4;

        // Complex refractive index
        dr::Complex<UnpolarizedSpectrum> base_eta_k(m_eta_base->eval(si, active),
                                                m_k_base->eval(si, active));
        if(dr::any_or<true>(selected_t)) {
            // ベース法線決定
            auto [a1, a2, a3, pa] = sampleRoute(wv1, m_pa, m_qa, m_ra, rng);
            auto [b1, b2, b3, pb] = sampleRoute(wv1, m_pb, m_qb, m_rb, rng);
            r1 = rand(rng);
            Mask elemA = (r1 <= 0.5f) && active;
            Mask elemB = (r1 >  0.5f) && active;
            n1 = dr::select(elemA, a1, b1);
            n2 = dr::select(elemA, a2, b2);
            n3 = dr::select(elemA, a3, b3);
            pathProb = dr::select(elemA, pa, pb);

            // n1反射
            r2 = rand2(rng);
            std::tie(m1, D1) = sample_distr_internal.sample(rotateVector(-wv1, n1, n), r2);
            m1 = rotateVector(m1, n, n1);
            F1 = fresnel_conductor(UnpolarizedSpectrum(dot(-wv1, m1)), base_eta_k);
            wv2 = reflect(-wv1, m1);
            // printf("%f, %f, %f, %f, %f\n", mag(wv1), mag(rotateVector(-wv1, n1, n)), mag(std::get<0>(sample_distr_internal.sample(rotateVector(-wv1, n1, n), r2))), mag(m1), mag(wv2));

            // n2反射
            r2 = rand2(rng);
            std::tie(m2, D2) = sample_distr_internal.sample(rotateVector(-wv2, n2, n), r2);
            m2 = rotateVector(m2, n, n2);
            F2 = fresnel_conductor(UnpolarizedSpectrum(dot(-wv2, m2)), base_eta_k);
            wv3 = reflect(-wv2, m2);

            // n3反射
            r2 = rand2(rng);
            std::tie(m3, D3) = sample_distr_internal.sample(rotateVector(-wv3, n3, n), r2);
            m3 = rotateVector(m3, n, n3);
            F3 = fresnel_conductor(UnpolarizedSpectrum(dot(-wv3, m3)), base_eta_k);
            wv4 = reflect(-wv3, m3);

            // n透過
            Float eta_v4_o, eta_o_v4, cos_theta_o;
            r2 = rand2(rng);
            std::tie(mo, Do) = sample_distr_surface.sample(wv4, r2);
            std::tie(Fo, cos_theta_o, eta_v4_o, eta_o_v4) = fresnel(dot(-wv4, mo), Float(m_eta));
            wo_rr = refract(-wv4, mo, cos_theta_o, eta_o_v4);
            // printf("(%f, %f, %f)check: rotate from(%f, %f, %f)to(%f, %f, %f) by <%f, %f, %f>\n", m2.x(), m2.y(), m2.z(), -wv2.x(), -wv2.y(), -wv2.z(), rotateVector(-wv2, n2, n).x(), rotateVector(-wv2, n2, n).y(), rotateVector(-wv2, n2, n).z(), n2.x(), n2.y(), n2.z());
            // printf("(%f, %f, %f, %f, %f, %f)  [%f, %f, %f, %f, %f]\n", 
            //         dot(si.wi,si.wi), dot(wv1, wv1), dot(wv2, wv2), dot(wv3, wv3), dot(wv4, wv4), dot(wo_rr, wo_rr),
            //         dot(mi, mi), dot(m1, m1), dot(m2, m2), dot(m3, m3), dot(mo, mo));
            // printf("[%f, %f, %f, %f, %f]\n",
            //     dot(mi, n), dot(m1, n1), dot(m2, n2), dot(m3, n3), dot(mo, n));
            // printf("[%f, %f, %f, %f, %f]  <%f, %f, %f> <%f, %f, %f>\n",
            //     dot(mi, n), dot(m1, n1), dot(m2, n2), dot(m3, n3), dot(mo, n),
            //     n2.x(), n2.y(), n2.z(), m2.x(), m2.y(), m2.z());
        }
        // incorrectPath: 経路決定後、次の面には裏面からしか入射しない場合
        Mask incorrectRRpath = !isfinite(Fo);
        F1 = dr::select(incorrectRRpath, 0.f, F1);
        F2 = dr::select(incorrectRRpath, 0.f, F2);
        F3 = dr::select(incorrectRRpath, 0.f, F3);
        Fo = dr::select(incorrectRRpath, 0.f, Fo);

        r1 = rand(rng);

        // ERA, 再帰反射可能な領域。拡散と再帰反射の割合を制御する
        Float era_cos_i = Frame3f::cos_theta(-wv1);
        auto [era_sin_phi, era_cos_phi] = Frame3f::sincos_phi(-wv1);
        Float era = ERA(era_cos_i, era_sin_phi, era_cos_phi);

        selected_rr = !selected_r && (r1<=(1.f-Fi)*(1.f-Fo) * era) && !incorrectRRpath && active;
        selected_d = !selected_r && !selected_rr && active;

        // printf("[%f, %f, %f]", r1, Fi, (1.f-Fi)*F1*F2*F3*(1.f-Fo));
        // printf("[%f, %f, %f, %f, %f](%f, %f)\n", Fi, F1, F2, F3, Fo, dot(-wv3, m3), dot(-wv4, mo));
        // printV(-wv3);printN(m3);printf(":%f\n", F3);
        /*
        printf("(%f, %f, %f)-s>(%f, %f, %f)-1>(%f, %f, %f)-2>(%f, %f, %f)-3>(%f, %f, %f)-s>(%f, %f, %f)  <%f, %f, %f> <%f, %f, %f> <%f, %f, %f>\n", 
                si.wi.x(), si.wi.y(), si.wi.z(),
                wv1.x(), wv1.y(), wv1.z(),
                wv2.x(), wv2.y(), wv2.z(),
                wv3.x(), wv3.y(), wv3.z(),
                wv4.x(), wv4.y(), wv4.z(),
                wo_rr.x(), wo_rr.y(), wo_rr.z(),
                m1.x(), m1.y(), m1.z(),
                m2.x(), m2.y(), m2.z(),
                m3.x(), m3.y(), m3.z());
        */
       
        // printN(mi);
        // printN(m1);
        // printN(m2);
        // printN(m3);
        // printN(mo);printf("\n");
        Float cos_theta_o = Frame3f::cos_theta(wo_rr);

        // それぞれの表面のBSDFSampleを作る
        // 表面
        if (dr::any_or<true>(selected_r)) {
            // printf("r");
            bs.wo[selected_r] = reflect(si.wi, mi);
            bs.pdf = dr::select(selected_r, Fi*Di*dr::rcp(4*dot(si.wi,mi)), bs.pdf);
            bs.eta = 1.f;
            bs.sampled_type = dr::select(selected_r, UInt32(+BSDFFlags::GlossyReflection), bs.sampled_type);
            bs.sampled_component = dr::select(selected_r, UInt32(0), bs.sampled_component);
            UnpolarizedSpectrum weight_r = 0.f;
            if (likely(m_sample_visible)) {
                weight_r = distr_surface.smith_g1(bs.wo, mi);
            } else {
                weight_r = distr_surface.G(si.wi, bs.wo, mi) * dot(si.wi, mi) /
                        (cos_theta_i * Frame3f::cos_theta(mi));
            }
            weight = dr::select(selected_r, weight_r, weight);
            weight[selected_r] *= m_surface_reflectance->eval(si, selected_r);
        }

        // 再帰
        if (dr::any_or<true>(selected_rr)) {
            // printf("R");
            // printV(si.wi);printV(wo_rr); printf("\n");
            Vector3f r3 = rand3c(rng);
            si.p[selected_rr] += Vector3f(randX*r3.x(), randY*r3.y(), randZ*r3.z());
            bs.wo[selected_rr] = wo_rr;
            bs.eta = dr::select(selected_rr, 1.f, bs.eta);
            bs.sampled_type = dr::select(selected_rr, UInt32(+BSDFFlags::GlossyTransmission | +BSDFFlags::GlossyReflection), bs.sampled_type);
            bs.sampled_component = dr::select(selected_rr, UInt32(1), bs.sampled_component);
            Spectrum weight_rr = 0.f;
            // ヤコビアン(bs.pdf)
            Float dwmi_dwi = dr::sqr(m_eta_air)*abs(dot(si.wi, mi)) / dr::sqr(m_eta_air*dot(si.wi, mi) + m_eta_mat*dot(wv1, mi));
            Float dwm1_dwv1 = dr::rcp(4.f * abs(dot(-wv1, m1)));
            Float dwm2_dwv2 = dr::rcp(4.f * abs(dot(-wv2, m2)));
            Float dwm3_dwv3 = dr::rcp(4.f * abs(dot(-wv3, m3)));
            Float dwmo_dwv4 = dr::sqr(m_eta_air)*abs(dot(wo_rr, mo)) / dr::sqr(m_eta_mat*dot(-wv4, mo) + m_eta_air*dot(wo_rr, mo));
            // bs.pdf = dr::select(selected_rr, 1.f, bs.pdf);
            bs.pdf = dr::select(selected_rr, (1.f-Fi)*(1.f-Fo) * dr::clamp(Di*D1*D2*D3*Do, 0.f, 1000000000000.f) * dwmi_dwi*dwm1_dwv1*dwm2_dwv2*dwm3_dwv3*dwmo_dwv4 * pathProb * era, bs.pdf);
            // bs.pdf = dr::select(selected_rr, (1.f-Fi)*F1*F2*F3*(1.f-Fo) *  dwmi_dwi*dwm1_dwv1*dwm2_dwv2*dwm3_dwv3*dwmo_dwv4 * Di*D1*D2*D3*Do * pathProb, bs.pdf);
            // bs.pdf = dr::select(selected_rr, (1.f-Fi)*F1*F2*F3*(1.f-Fo) * Di*D1*D2*D3*Do * dwmi_dwi*dwm1_dwv1*dwm2_dwv2*dwm3_dwv3*dwmo_dwv4 * pathProb, bs.pdf);
            // bs.pdf = dr::select(selected_rr, 1.f* pathProb, bs.pdf);
            // printf("<F:%f, D:%f, J:%f, p:%f, pdf:%f>\n",  (1.f-Fi)*F1*F2*F3*(1.f-Fo), min(Di*D1*D2*D3*Do, 1000000000), dwmi_dwi*dwm1_dwv1*dwm2_dwv2*dwm3_dwv3*dwmo_dwv4, pathProb, bs.pdf);
            // printV(si.wi);printV(wv1);printV(wv2);printV(wv3);printV(wv4);printV(wo_rr);
            // printf("<%f, %f, %f, %f, %f, -> %f>\n",  Di, D1, D2, D3, Do, Di*D1*D2*D3*Do);
            //////////////////////////// pdf /////////////////////////////////
            if (likely(m_sample_visible)) {
                weight_rr = distr_surface.smith_g1(wv1, mi) *
                            distr_internal.smith_g1(rotateVector(wv2, n1, n), rotateNormal(m1, n1, n)) *
                            distr_internal.smith_g1(rotateVector(wv3, n2, n), rotateNormal(m2, n2, n)) *
                            distr_internal.smith_g1(rotateVector(wv4, n3, n), rotateNormal(m3, n3, n)) *
                            distr_surface.smith_g1(wo_rr, mo);
                // printf("<%f, %f>\t", distr_surface.smith_g1(wv1, mi) *
                //             distr_internal.smith_g1(rotateVector(wv2, n1, n), rotateNormal(m1, n1, n)) *
                //             distr_internal.smith_g1(rotateVector(wv3, n2, n), rotateNormal(m2, n2, n)) *
                //             distr_internal.smith_g1(rotateVector(wv4, n3, n), rotateNormal(m3, n3, n)) *
                //             distr_surface.smith_g1(wo_rr, mo), weight_rr);
            } else {
                weight_rr = distr_surface.G(si.wi, wv1, mi) *
                            distr_internal.G(rotateVector(-wv1, n1, n), rotateVector(wv2, n1, n), rotateNormal(m1, n1, n)) *
                            distr_internal.G(rotateVector(-wv2, n2, n), rotateVector(wv3, n2, n), rotateNormal(m2, n2, n)) *
                            distr_internal.G(rotateVector(-wv3, n3, n), rotateVector(wv4, n3, n), rotateNormal(m3, n3, n)) *
                            distr_surface.G(-wv4, wo_rr, mo) *
                            abs(dot(si.wi,mi)*dot(wv1,m1)*dot(wv2,m2)*dot(wv3,m3)*dot(wv4,mo)) *
                            dr::rcp(cos_theta_i*cos_theta_v1 * dot(wv1,n1)*dot(wv2,n1) * dot(wv2,n2)*dot(wv3,n2) * dot(wv3,n3)*dot(wv4,n3) * dot(wv4,n)*cos_theta_o);
            }
            weight = dr::select(selected_rr, weight_rr, weight);
            // printf("[%f, %f, %f, %f, %f -> %f], %f\n", distr_surface.smith_g1(wv1, mi),
            //                 distr_internal.smith_g1(rotateVector(wv2, n1, n), rotateNormal(m1, n1, n)),
            //                 distr_internal.smith_g1(rotateVector(wv3, n2, n), rotateNormal(m2, n2, n)),
            //                 distr_internal.smith_g1(rotateVector(wv4, n3, n), rotateNormal(m3, n3, n)),
            //                 distr_surface.smith_g1(wo_rr, mo), distr_surface.smith_g1(wv1, mi) *
            //                 distr_internal.smith_g1(rotateVector(wv2, n1, n), rotateNormal(m1, n1, n)) *
            //                 distr_internal.smith_g1(rotateVector(wv3, n2, n), rotateNormal(m2, n2, n)) *
            //                 distr_internal.smith_g1(rotateVector(wv4, n3, n), rotateNormal(m3, n3, n)) *
            //                 distr_surface.smith_g1(wo_rr, mo), weight_rr);
            // weight = dr::select(selected_rr, 1.f, weight);
            weight[selected_rr] *= m_surface_reflectance->eval(si, selected_rr)*m_surface_reflectance->eval(si, selected_rr);
            weight[selected_rr] *= m_internal_reflectance->eval(si, selected_rr)*m_internal_reflectance->eval(si, selected_rr)*m_internal_reflectance->eval(si, selected_rr);
            
        }

        // 拡散
        if (dr::any_or<true>(selected_d)) {
            // printf("D");
            r2 = rand2(rng);
            Vector3f wo_diffuse = warp::square_to_cosine_hemisphere(r2);
            bs.wo[selected_d] = wo_diffuse;
            bs.pdf = dr::select(selected_d, (1.f-Fi) * (era * Fo * pathProb + (1.f-era)) * warp::square_to_cosine_hemisphere_pdf(wo_diffuse), bs.pdf);
            bs.eta = dr::select(selected_d, 1.f, bs.eta);
            bs.sampled_type = dr::select(selected_d, UInt32(+BSDFFlags::DiffuseReflection), bs.sampled_type);
            bs.sampled_component = dr::select(selected_d, UInt32(2), bs.sampled_component);
            weight = dr::select(selected_d, 1.f, weight);
        }
        // weight = 0.1f;
        // printf("<%f, %f>\n",(1.f-Fi)*F1*F2*F3*(1.f-Fo), (1.f-Fi)*(1.f-F1*F2*F3*(1.f-Fo)));
        
        return {bs, weight};
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);
        MicrofacetDistribution distr_surface(m_type,
                                             m_alpha_u_surface->eval_1(si, active),
                                             m_alpha_v_surface->eval_1(si, active),
                                             m_sample_visible);
        // return 1.f;
        // BRDF, activeチェック
        Spectrum value = 0.f;
        bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        bool has_retroreflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        has_retroreflection &= ctx.is_enabled(BSDFFlags::GlossyTransmission, 1);
        bool has_diffusereflection = ctx.is_enabled(BSDFFlags::DiffuseReflection, 2);
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float cos_theta_o = Frame3f::cos_theta(wo);
        active &= (cos_theta_i>0.f) && (cos_theta_o>0.f);
        if (unlikely(dr::none_or<false>(active)))
            return 0.f;
        // Complex refractive index
        dr::Complex<UnpolarizedSpectrum> base_eta_k(m_eta_base->eval(si, active),
                                                m_k_base->eval(si, active));
        // 12経路の内部相互作用の評価
        Path p[12] = {
            Path(n, m_pa, m_qa, m_ra, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pa, m_ra, m_qa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qa, m_pa, m_ra, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qa, m_ra, m_pa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_ra, m_pa, m_qa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_ra, m_qa, m_pa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pb, m_qb, m_rb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pb, m_rb, m_qb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qb, m_pb, m_rb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qb, m_rb, m_pb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_rb, m_pb, m_qb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_rb, m_qb, m_pb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
        };

        Spectrum brdf_r  = 0.f;
        Spectrum brdf_rr = 0.f;
        Spectrum brdf_d  = 0.f;
        // 表面
        Vector3f m = normalize(si.wi + wo);
        Float F = std::get<0>(fresnel(dot(si.wi, m), Float(m_eta)));
        Float G = distr_surface.G(si.wi, wo, m);
        Float D = distr_surface.eval(m);
        brdf_r = F*G*D* 0.25f*dr::rcp(cos_theta_i); // cos_theta_o打ち消し
        // printf("(%f, %f, %f, %f)\t", F, G, D, brdf_r);

        // 再帰・拡散
        for (int i=0; i<12; i++) {
            auto [v_rr, v_d] = p[i].eval(si, m_alpha_u_surface, m_alpha_v_surface, m_alpha_u_internal, m_alpha_v_internal, m_type, active, m_sample_visible);
            brdf_rr += v_rr;
            brdf_d += v_d;
        }
        brdf_rr *= cos_theta_o;
        brdf_d *= diffuseFactor * cos_theta_o;

        brdf_r *= m_surface_reflectance->eval(si, active);
        brdf_rr *= m_surface_reflectance->eval(si, active)*m_surface_reflectance->eval(si, active);
        brdf_rr *= m_internal_reflectance->eval(si, active)*m_internal_reflectance->eval(si, active)*m_internal_reflectance->eval(si, active);

        // 足して返す
        value = dr::select(active, brdf_d + brdf_rr + brdf_r, value);
        // printf("(%f, %f, %f)\n", brdf_r, brdf_rr, brdf_d);
        // Mask invalid = !isfinite(brdf_d) || !isfinite(brdf_rr) || !isfinite(brdf_d);
        // if (dr::any_or<true>(invalid)) { printf("invalid:(%f, %f, %f)\t", brdf_d, brdf_rr, brdf_r);}
        return (active, value, 0.f);
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);
        // return 1.f;
        // activeチェック
        bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        bool has_retroreflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        has_retroreflection &= ctx.is_enabled(BSDFFlags::GlossyTransmission, 1);
        bool has_diffusereflection = ctx.is_enabled(BSDFFlags::DiffuseReflection, 2);
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float cos_theta_o = Frame3f::cos_theta(wo);
        active &= (cos_theta_i>0.f) && (cos_theta_o>0.f);
        
        if (unlikely(dr::none_or<false>(active)))
            return 0.f;
        // Distributions, Walter's Trick
        MicrofacetDistribution distr_surface(   m_type,
                                                m_alpha_u_surface->eval_1(si, active),
                                                m_alpha_v_surface->eval_1(si, active),
                                                m_sample_visible);
        // Walter's trick
        if (unlikely(!m_sample_visible))
            distr_surface.scale_alpha(1.2f - .2f * sqrt(abs(Frame3f::cos_theta(si.wi))));
        
        Float prob_r  = 0.f;
        Float prob_rr = 0.f;
        Float prob_d  = 0.f;
        // 表面
        Vector3f m = normalize(si.wi + wo);
        prob_r = distr_surface.pdf(si.wi, m) * std::get<0>(fresnel(dot(si.wi, m), Float(m_eta))) * 0.25f *dr::rcp(abs_dot(si.wi, m));

        // 再帰・拡散
        // Complex refractive index
        dr::Complex<UnpolarizedSpectrum> base_eta_k(m_eta_base->eval(si, active),
                                                m_k_base->eval(si, active));
        Path p[12] = {
            Path(n, m_pa, m_qa, m_ra, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pa, m_ra, m_qa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qa, m_pa, m_ra, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qa, m_ra, m_pa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_ra, m_pa, m_qa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_ra, m_qa, m_pa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pb, m_qb, m_rb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pb, m_rb, m_qb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qb, m_pb, m_rb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qb, m_rb, m_pb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_rb, m_pb, m_qb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_rb, m_qb, m_pb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
        };
        for (int i=0; i<12; i++) {
            auto [v_rr, v_d] = p[i].pdf(si, m_alpha_u_surface, m_alpha_v_surface, m_alpha_u_internal, m_alpha_v_internal, m_type, active, m_sample_visible);
            prob_rr += v_rr;
            prob_d += v_d;
        }
        // prob_d *= diffuseFactor * warp::square_to_cosine_hemisphere_pdf(wo);
        prob_d *= diffuseFactor;
        Mask invalid = !isfinite(prob_r) || !isfinite(prob_rr) || !isfinite(prob_d);
        // if (dr::any_or<true>(invalid)) { printf("invalid:[%f, %f, %f]\t", prob_d, prob_rr, prob_r);}
        // 足して返す
        // printf("[%f, %f, %f]\n", prob_r, prob_rr, prob_d);       
        return dr::select(active && !invalid, prob_r + prob_rr + prob_d, 0.f);

    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);
        MicrofacetDistribution distr_surface(m_type,
                                             m_alpha_u_surface->eval_1(si, active),
                                             m_alpha_v_surface->eval_1(si, active),
                                             m_sample_visible);
        bool has_reflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 0);
        bool has_retroreflection = ctx.is_enabled(BSDFFlags::GlossyReflection, 1);
        has_retroreflection &= ctx.is_enabled(BSDFFlags::GlossyTransmission, 1);
        bool has_diffusereflection = ctx.is_enabled(BSDFFlags::DiffuseReflection, 2);
        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float cos_theta_o = Frame3f::cos_theta(wo);
        active &= (cos_theta_i>0.f) && (cos_theta_o>0.f);
        if (unlikely(dr::none_or<false>(active)))
            return {0.f, 0.f};
        
        Spectrum brdf_r  = 0.f;
        Spectrum brdf_rr = 0.f;
        Spectrum brdf_d  = 0.f;
        Float prob_r  = 0.f;
        Float prob_rr = 0.f;
        Float prob_d  = 0.f;
        Vector3f m = normalize(si.wi + wo);

        // surface
        Float F = std::get<0>(fresnel(dot(si.wi, m), Float(m_eta)));
        Float G = distr_surface.G(si.wi, wo, m);
        Float D = distr_surface.eval(m);
        brdf_r = F*G*D* 0.25f*dr::rcp(cos_theta_i); // cos_theta_o打ち消し
        prob_r = distr_surface.pdf(si.wi, m) * std::get<0>(fresnel(dot(si.wi, m), Float(m_eta))) * 0.25f *dr::rcp(abs_dot(si.wi, m));

        //retroreflect / diffuse
        Path p[12] = {
            Path(n, m_pa, m_qa, m_ra, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pa, m_ra, m_qa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qa, m_pa, m_ra, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qa, m_ra, m_pa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_ra, m_pa, m_qa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_ra, m_qa, m_pa, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pb, m_qb, m_rb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_pb, m_rb, m_qb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qb, m_pb, m_rb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_qb, m_rb, m_pb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_rb, m_pb, m_qb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
            Path(n, m_rb, m_qb, m_pb, si.wi, wo, m_eta_air, m_eta_mat, m_eta_base, m_k_base),
        };
        for (int i=0; i<12; i++) {
            auto [v_rr_brdf, v_d_brdf] = p[i].eval(si, m_alpha_u_surface, m_alpha_v_surface, m_alpha_u_internal, m_alpha_v_internal, m_type, active, m_sample_visible);
            auto [v_rr_prob, v_d_prob] = p[i].pdf( si, m_alpha_u_surface, m_alpha_v_surface, m_alpha_u_internal, m_alpha_v_internal, m_type, active, m_sample_visible);
            brdf_rr += v_rr_brdf;
            brdf_d += v_d_brdf;
            prob_rr += v_rr_prob;
            prob_d += v_d_prob;
        }

        brdf_rr *= cos_theta_o;
        brdf_d *= diffuseFactor * cos_theta_o;
        prob_d *= diffuseFactor;

        brdf_r *= m_surface_reflectance->eval(si, active);
        brdf_rr *= m_surface_reflectance->eval(si, active)*m_surface_reflectance->eval(si, active);
        brdf_rr *= m_internal_reflectance->eval(si, active)*m_internal_reflectance->eval(si, active)*m_internal_reflectance->eval(si, active);
        Mask invalid = !isfinite(prob_r) || !isfinite(prob_rr) || !isfinite(prob_d);

        auto brdf = dr::select(active, brdf_d + brdf_rr + brdf_r, 0.f);
        auto prob = dr::select(active && invalid, prob_r + prob_rr + prob_d, 0.f);

        return {brdf, prob};

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
    ref<Texture> m_surface_reflectance;
    ref<Texture> m_internal_reflectance;
    Float randX, randY, randZ;
};

MI_IMPLEMENT_CLASS_VARIANT(RoughRetroreflector, BSDF)
MI_EXPORT_PLUGIN(RoughRetroreflector, "Rough retroreflector")
NAMESPACE_END(mitsuba)
