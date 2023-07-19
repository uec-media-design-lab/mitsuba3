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
        m_pa = Normal3f(2*math::InvSqrtSix<Float>, 0, math::InvSqrtThree<Float>);                         // →
        m_qa = Normal3f(-math::InvSqrtSix<Float>, math::InvSqrtTwo<Float>, math::InvSqrtThree<Float>);    // ↖
        m_ra = Normal3f(-math::InvSqrtSix<Float>, -math::InvSqrtTwo<Float>, math::InvSqrtThree<Float>);   // ↙
        m_pb = Normal3f(-2*math::InvSqrtSix<Float>, 0, math::InvSqrtThree<Float>);                       // ←
        m_qb = Normal3f(math::InvSqrtSix<Float>, math::InvSqrtTwo<Float>, math::InvSqrtThree<Float>);    // ↗
        m_rb = Normal3f(math::InvSqrtSix<Float>, -math::InvSqrtTwo<Float>, math::InvSqrtThree<Float>);   // ↘

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
    
        // sample_visible
        m_sample_visible = props.bool_("sample_visible", true);

        // diffuseFactor
        Float F_d = 0.f;
        if (m_eta<1.f) {
            F_d = 0.919317 - 3.4793*m_eta + 6.75335*m_eta*m_eta - 7.80989*m_eta*m_eta*m_eta + 4.98554*m_eta*m_eta*m_eta*m_eta - 1.36881*m_eta*m_eta*m_eta*m_eta*m_eta;
        } else {
            F_d = -9.23372 + 22.2272*m_eta - 20.9292*m_eta*m_eta + 10.2291*m_eta*m_eta*m_eta - 2.54396*m_eta*m_eta*m_eta*m_eta + 0.254913*m_eta*m_eta*m_eta*m_eta*m_eta;
        }
        diffuseFactor = math::InvPi<Float> * sqr(m_inv_eta) / (1.f - F_d);
        
        // surface normal(ideal)
        n = Normal3f(0.f, 0.f, 1.f);

        // reflectance
        m_surface_reflectance   = props.float_("surface_reflectance", 1.f);
    
        m_internal_reflectance   = props.float_("internal_reflectance", 1.f);
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

    
    // print vector(or normal)(s)
    void printV(const Vector3f v) const {
        printf("(%f, %f, %f)\t", v.x(), v.y(), v.z());
    }
    void printN(const Normal3f n) const {
        printf("<%f, %f, %f>\t", n.x(), n.y(), n.z());
    }

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
        Float inner_a = clamp(dot(-wi, a), 0.f, 1.f);
        Float inner_b = clamp(dot(-wi, b), 0.f, 1.f);
        Float inner_c = clamp(dot(-wi, c), 0.f, 1.f);
        Float invDenominator = 1.f / (inner_a + inner_b + inner_c);
        n1 =    select(Mask(r2.x() <= inner_a*invDenominator), a,
                select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), b, c));
        prob =  select(Mask(r2.x() <= inner_a*invDenominator), inner_a,
                select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), inner_b, inner_c)) * invDenominator * 0.5f;

        // 2nd, 3rd surfaces
        nt1 =   select(Mask(r2.x() <= inner_a*invDenominator), b, 
                select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), a, a));
        nt2 =   select(Mask(r2.x() <= inner_a*invDenominator), c, 
                select(Mask(r2.x() <= (inner_a+inner_b)*invDenominator), c, b));
        vt = reflect(-wi, n1);
        Float inner_nt1 = clamp(dot(-vt, nt1), 0.f, 1.f);
        Float inner_nt2 = clamp(dot(-vt, nt2), 0.f, 1.f);
        invDenominator = 1.f / (inner_nt1 + inner_nt2);

        n2 =            select(Mask(r2.y() <= inner_nt1*invDenominator), nt1, nt2);
        n3 =            select(Mask(r2.y() <= inner_nt1*invDenominator), nt2, nt1);
        prob = prob *   select(Mask(r2.y() <= inner_nt1*invDenominator), inner_nt1, inner_nt2) * invDenominator;

        return std::forward_as_tuple(n1, n2, n3, prob);
    }

    // ERA
    Float ERA(const Float cos_theta_i, const Float sin_phi_i, const Float cos_phi_i) const {
        Float x = rad_to_deg(acos(abs(cos_theta_i)));
        Float x_phi = mulsign(rad_to_deg(acos(cos_phi_i)), sin_phi_i) + 360.f;  // 0 -- 360 [deg]
        Float x_phi_blend = abs(fmod(x_phi, 60.f) - 30.f) * 0.1f * 0.33333333333333333333f;
        Float val_upper = ERA_upperBoundary(x);
        Float val_lower = ERA_lowerBoundary(x);
        // printf("phi:%f\tphi:%f\tphimod:%f\tphimodabs30:%f\tblend:%f\n",rad_to_deg(acos(cos_phi_i)), x_phi, fmod(x_phi, 60.f), abs(fmod(x_phi, 60.f) - 30.f), x_phi_blend);
        // printf("cos:%f\tsin:%f\tphi:%f\tphi':%f\tphi-blend:%f\n", cos_phi_i, sin_phi_i, x_phi, rad_to_deg(acos(cos_phi_i)), x_phi_blend);
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
        return select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
    }
    Float ERA_upperBoundary(const Float x) const {
        Float a = -0.353394;
        Float b = 0.055244;
        Float c = -1.632606;
        Float d = 0.322462;
        Float val = a*atan(b*x +c)+d;
        // printf("x:%f, val:%f\n", x, val);
        return select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
    }

    // class
    class Path {
        // 
        public:
            // MTS_IMPORT_BASE(BSDF, m_flags, m_components)
            // MTS_IMPORT_TYPES(Texture, MicrofacetDistribution)
            MTS_IMPORT_TYPES()
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
                Float inner_1st_n1 = clamp(dot(-muv1, m1), 0.f, 1.f);
                Float inner_1st_n2 = clamp(dot(-muv1, m2), 0.f, 1.f);
                Float inner_1st_n3 = clamp(dot(-muv1, m3), 0.f, 1.f);
                Float invDenominator1st = 1.f / (inner_1st_n1 + inner_1st_n2 + inner_1st_n3);

                prob = prob * inner_1st_n1*invDenominator1st;

                Float inner_2nd_n2 = clamp(dot(-muv2, m2), 0.f, 1.f);
                Float inner_2nd_n3 = clamp(dot(-muv2, m3), 0.f, 1.f);
                Float invDenominator2nd = 1.f / (inner_2nd_n2 + inner_2nd_n3);
                
                prob = prob * inner_2nd_n2*invDenominator2nd;
                
                valid = Mask(inner_1st_n1 + inner_1st_n2 + inner_1st_n3 > 0.f)&& valid;
                valid = Mask(inner_2nd_n2 + inner_2nd_n3>0.f)&& valid;
                valid = Mask(isfinite(prob)) && valid;
                return select(valid, prob, 0.f);
            }

            // class - ERA
            Float ERA(const Float cos_theta_i, const Float sin_phi_i, const Float cos_phi_i) const {
            Float x = rad_to_deg(acos(abs(cos_theta_i)));
            Float x_phi = mulsign(rad_to_deg(acos(cos_phi_i)), sin_phi_i) + 360.f;  // 0 -- 360 [deg]
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
            return select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
        }
        Float ERA_upperBoundary(const Float x) const {
            Float a = -0.353394;
            Float b = 0.055244;
            Float c = -1.632606;
            Float d = 0.322462;
            Float val = a*atan(b*x +c)+d;
            // printf("x:%f, val:%f\n", x, val);
            return select(Mask(val>=0), val, 0.f);  // Max(val, 0.f)
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
            Complex<UnpolarizedSpectrum> base_eta_k(eta_base->eval(si, active),
                                                    k_base->eval(si, active));
            Float Fi = std::get<0>(fresnel(dot(-mui, mi), Float(reta)));
            auto F1 = fresnel_conductor(UnpolarizedSpectrum(dot(-muv1, m1)), base_eta_k);
            auto F2 = fresnel_conductor(UnpolarizedSpectrum(dot(-muv2, m2)), base_eta_k);
            auto F3 = fresnel_conductor(UnpolarizedSpectrum(dot(-muv3, m3)), base_eta_k);
            Float Fo = std::get<0>(fresnel(dot(muo, mo), Float(reta)));
            // incorrect Path check
            // printf("(%f:  %f, %f, %f, %f, %f)\n", pathProb, Fi, F1, F2, F3, Fo);
            Mask incorrectRRpath = !isfinite(Fo);
            Fi = select(incorrectRRpath, 0.f, Fi);
            F1 = select(incorrectRRpath, 0.f, F1);
            F2 = select(incorrectRRpath, 0.f, F2);
            F3 = select(incorrectRRpath, 0.f, F3);
            Fo = select(incorrectRRpath, 0.f, Fo);

            // G
            Float Gi = distr_surface.G(-mui, muv1, mi);
            Float G1 = distr_internal.G(rotateVector(-muv1,n1,ni), rotateVector(muv2,n1,ni), rotateNormal(m1,n1,ni));
            Float G2 = distr_internal.G(rotateVector(-muv2,n2,ni), rotateVector(muv3,n2,ni), rotateNormal(m2,n2,ni));
            Float G3 = distr_internal.G(rotateVector(-muv3,n3,ni), rotateVector(muv4,n3,ni), rotateNormal(m3,n3,ni));
            Float Go = distr_surface.G(-muv4, muo, mo);

            // J
            Float Ji = abs(ea*ea*dot(mui, mi))*abs(dot(muv1, mi)) / (abs(dot(mui, ni))*abs(dot(muv1, ni)) * sqr(ea*dot(-mui,mi) + em*dot(muv1,mi)));
            Float J1 = 0.25f / (abs(dot(muv1,n1)*dot(muv2,n1)));
            Float J2 = 0.25f / (abs(dot(muv2,n2)*dot(muv3,n2)));
            Float J3 = 0.25f / (abs(dot(muv3,n3)*dot(muv4,n3)));
            Float Jo = abs(em*em*dot(muv4, mo))*abs(dot(muo, mo)) / (abs(dot(muv4, no))*abs(dot(muo, no)) * sqr(em*dot(-muv4,mo) + ea*dot(muo,mo)));

            // ERA
            Float era_cos_i = Frame3f::cos_theta(-muv1);
            auto [era_sin_phi, era_cos_phi] = Frame3f::sincos_phi(-muv1);
            Float era = ERA(era_cos_i, era_sin_phi, era_cos_phi);

            // RetroReflection
            // val_rr = pathProb * (1.f-Fi)*F1*F2*F3*(1.f-Fo) * pow(Gi*G1*G2*G3*Go * Di*D1*D2*D3*Do * Ji*J1*J2*J3*Jo, 0.2);
            val_rr = select(!incorrectRRpath, 0.f, pathProb * (1.f-Fi)*F1*F2*F3*(1.f-Fo) * Gi*G1*G2*G3*Go * Di*D1*D2*D3*Do * Ji*J1*J2*J3*Jo * era);
            // printf("[P:%f, F:%f, G:%f, D:%f, J:%f]\n", pathProb, (1.f-Fi)*F1*F2*F3*(1.f-Fo), Gi*G1*G2*G3*Go, pow(Di*D1*D2*D3*Do, 0.2), Ji*J1*J2*J3*Jo);

            // Diffuse
            // diffuseFactorを除いた部分
            val_d = select(isfinite(pathProb), (1.f-Fi)* (era * pathProb * (1.f-F1*F2*F3*(1.f-Fo)) + (1.f-era)), 0.f);
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
            Complex<UnpolarizedSpectrum> base_eta_k(eta_base->eval(si, active),
                                                    k_base->eval(si, active));
            Float Fi = std::get<0>(fresnel(dot(-mui, mi), Float(reta)));
            Float F1 = 1.f;
            Float F2 = 1.f;
            Float F3 = 1.f;
            Float Fo = std::get<0>(fresnel(dot(muo, mo), Float(reta)));
            // incorrect Path check
            Mask incorrectRRpath = !isfinite(Fo);
            F1 = select(incorrectRRpath, 0.f, F1);
            F2 = select(incorrectRRpath, 0.f, F2);
            F3 = select(incorrectRRpath, 0.f, F3);
            Fo = select(incorrectRRpath, 0.f, Fo);

            // J
            Float Ji = em*em*abs(dot(muv1, mi)) / sqr(ea*dot(-mui,mi) + em*dot(muv1,mi));
            Float J1 = 0.25f / abs(dot(muv2,n1));
            Float J2 = 0.25f / abs(dot(muv3,n2));
            Float J3 = 0.25f / abs(dot(muv4,n3));
            Float Jo = ea*ea*abs(dot(muo, mo)) / sqr(em*dot(-muv4,mo) + ea*dot(muo,mo));

            // ERA
            Float era_cos_i = Frame3f::cos_theta(-muv1);
            auto [era_sin_phi, era_cos_phi] = Frame3f::sincos_phi(-muv1);
            Float era = ERA(era_cos_i, era_sin_phi, era_cos_phi);

            // RetroReflection
            prob_rr = select(!incorrectRRpath, 0.f, pathProb * (1.f-Fi)*F1*F2*F3*(1.f-Fo) * Di*D1*D2*D3*Do * Ji*J1*J2*J3*Jo * era);

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
