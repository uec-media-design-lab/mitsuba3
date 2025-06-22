#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/texture.h>

#include <drjit/dynamic.h>
#include <mitsuba/core/distr_2d.h>
#include <mitsuba/core/distr_1d.h>

#include <mitsuba/render/sampler.h>

#include <mitsuba/render/microfacet.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <random>
#include <iterator>
#include <string>
#include "nlohmann/json.hpp"
#include <vector>
#include <drjit/dynamic.h>
using namespace std;
using json = nlohmann::json;


NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class sinc_shift4 final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    using Index = dr::uint32_array_t<Float>;

    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    sinc_shift4(const Properties &props) : Base(props)
    {
        // ------ 表面反射光関係 -------
        // 表面の反射率
        if (props.has_property("specular_reflectance"))
            m_specular_reflectance   = props.texture<Texture>("specular_reflectance", 1.f);

        // 表面の透過率
        if (props.has_property("specular_transmittance"))
            m_specular_transmittance = props.texture<Texture>("specular_transmittance", 1.f);

        // 屈折率
        ScalarFloat int_ior = lookup_ior(props, "int_ior", "bk7");
        ScalarFloat ext_ior = lookup_ior(props, "ext_ior", "air");
        if (int_ior < 0.f || ext_ior < 0.f || int_ior == ext_ior)
            Throw("The interior and exterior indices of "
                  "refraction must be positive and differ!");
        m_eta = int_ior / ext_ior;
        m_inv_eta = ext_ior / int_ior;

        // Microfacet分布
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

        // 可視（Smith）
        m_sample_visible = props.get<bool>("sample_visible", true);

        // 表面の粗さ
        if (props.has_property("alpha_u") || props.has_property("alpha_v")) {
            if (!props.has_property("alpha_u") || !props.has_property("alpha_v"))
                Throw("Microfacet model: both 'alpha_u' and 'alpha_v' must be specified.");
            if (props.has_property("alpha"))
                Throw("Microfacet model: please specify"
                      "either 'alpha' or 'alpha_u'/'alpha_v'.");
            m_alpha_u = props.texture<Texture>("alpha_u");
            m_alpha_v = props.texture<Texture>("alpha_v");
        } else {
            m_alpha_u = m_alpha_v = props.texture<Texture>("alpha", 0.1f);
        }

        BSDFFlags extra = (m_alpha_u != m_alpha_v) ? BSDFFlags::Anisotropic : BSDFFlags(0);
        m_components.push_back(BSDFFlags::GlossyReflection | BSDFFlags::FrontSide |
                               BSDFFlags::BackSide | extra);
        m_components.push_back(BSDFFlags::GlossyTransmission | BSDFFlags::FrontSide |
                               BSDFFlags::BackSide | BSDFFlags::NonSymmetric | extra);
        m_flags = m_components[0] | m_components[1];

        parameters_changed();

        // ------ 再帰反射光関係 -------
        m_flags = m_flags | BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        // 反射率や拡がり角パラメータ、コーナーサイズなど
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        m_a = props.texture<Texture>("a", 1.f);
        m_angle = props.texture<Texture>("angle", 0.f);
        m_cornersize = props.texture<Texture>("cornersize", 0.f);
        m_shiftoffset = props.texture<Texture>("shiftoffset", 0.f);

        m_sample_visible = props.get<bool>("sample_visible", true);

        // sincLUTのサイズ（１つの角度に対して）
        M = 1; N = 30000;
        offset = 0.0;

        m_isEstimation = 1;

        // 光線シフト分布の作成
        make_rayshift_array();
    }

    void parameters_changed(const std::vector<std::string> &/*keys*/ = {}) override {
        m_inv_eta = dr::rcp(m_eta);
        dr::make_opaque(m_eta, m_inv_eta);
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

    void make_sincarray1(float a, int expo, int angle)
    {
        ScalarFloat data[N*10];
        ScalarVector2u size(N, 10);
        int idx = 0;
        float sums[10];

        for (int i = 0; i < 10; ++i)
        {
            a = alist[i];
            sums[i] = 0.0f;
            for (int j = 0; j < N; ++j) {
                ScalarFloat arg = (a*300/dr::cos(offset)) * dr::tan((dr::Pi<Float>*j/(2.0f*N) - offset));
                ScalarFloat sincvalue = dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo);

                sincvalue = dr::select(
                    dr::isfinite(sincvalue),
                    dr::abs(sincvalue),
                    ScalarFloat(1.0)
                );
                if (j >= 28000) sincvalue = 0.0f;
                data[idx++] = sincvalue;
                sums[i] += sincvalue;
            }
            printf("a[%d] = %f, sum = %f", i, a, sums[i]);
        }

        // 正規化（ごり押し）
        idx = 0;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < N; ++j)
            {
                data[idx++] /= sums[i];
                if (j < 20) printf("data[%d] = %f\n", i*N+j, data[i*N+j]);
            }
        }

        printf("size in make x = %d, y = %d\n", size.x(), size.y());
        DiscreteDistribution2D<Float, 2> dd(&data[0], size);
        m_pdfdata = dd;

        float reflist_float[10] = {0.3149740397930145, 0.31731319427490234, 0.303494393825531, 0.2993064522743225, 0.2851271331310272, 0.26412859559059143, 0.22004221379756927, 0.21191012859344482, 0.20523042976856232};
        ScalarFloat data2[10];
        for (int i = 0; i < 9; ++i)
            data2[i] = reflist_float[i];
        DiscreteDistribution<Float> rl(&data2[0], 10);
        reflist = rl;
    }

    // TODO 角度ごとのLUTを保存するように書き替える
    void make_rayshift_array() {
        vector<vector<float>> m_rayshift_list;
        std::string filename = "/home/sugawara.ryo/makeLUT/2025.03.17diff/LUTs/rayshift_LUT.json";
        ifstream ifs(filename.c_str());
        if (ifs.good())
        {
            json m_json;
            ifs >> m_json;

            for (const auto& json : m_json) {
                vector<float> d = json.get<vector<float>>();
                m_rayshift_list.push_back(d);
            }
        }
        else
        {
            cout << "ファイルの読み込みに失敗しました" << endl;
        }

        int num_i = 10;
        int num_j = m_rayshift_list[0].size();
        ScalarFloat data[num_i*num_j];
        
        int idx = 0;
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < num_j; ++j) {
                if (j == 0) printf("m_rayshift_list[%d][%d] = %f\n", i, j, m_rayshift_list[i][j]);
                data[idx++] = m_rayshift_list[i][j];
            }
        }
        printf("num_i = %d, num_j = %d, idx = %d\n", num_i, num_j, idx);
        ScalarVector2u size(num_j, num_i);
        DiscreteDistribution2D<Float, 2> dd(&data[0], size);
        m_rayshiftLUT_list = dd;
        rayshift_N = num_j;
    }

    void traverse(TraversalCallback *callback) override
    {
        // 表面反射光関係
        callback->put_parameter("eta", m_eta, ParamFlags::Differentiable | ParamFlags::Discontinuous);
        if (!has_flag(m_flags, BSDFFlags::Anisotropic))
            callback->put_object("alpha",                  m_alpha_u.get(),                ParamFlags::Differentiable | ParamFlags::Discontinuous);
        else {
            callback->put_object("alpha_u",                m_alpha_u.get(),                ParamFlags::Differentiable | ParamFlags::Discontinuous);
            callback->put_object("alpha_v",                m_alpha_v.get(),                ParamFlags::Differentiable | ParamFlags::Discontinuous);
        }
        if (m_specular_reflectance)
            callback->put_object("specular_reflectance",   m_specular_reflectance.get(),   +ParamFlags::Differentiable);
        if (m_specular_transmittance)
            callback->put_object("specular_transmittance", m_specular_transmittance.get(), +ParamFlags::Differentiable);
        
        // 再帰反射光関係
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("angle", m_angle.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("cornersize", m_cornersize.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("shiftoffset", m_shiftoffset.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);

        angle = static_cast<int>(m_angle.get()->max());

        angle_index = angle / 5.0;
        a = m_a.get()->max(); reflectance = m_reflectance.get()->max();
        
        make_sincarray1(a, 2, angle);
    }

    float generate_random() const
    {
        random_device seed_gen;
        static std::mt19937 generator(seed_gen());
        static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
        return distribution(generator);
    }


    mitsuba::PCG32<UInt32> setRandomGenerator(const Float seed) const {
        mitsuba::PCG32<UInt32> rng(1, PCG32_DEFAULT_STATE, seed);
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


    std::pair<BSDFSample3f, Spectrum> sample(
        const BSDFContext &ctx,
        const SurfaceInteraction3f &si,
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
        Float theta_i = dr::acos(cos_theta_i);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};

        auto rng = setRandomGenerator(sample1*10000000);

        // ----------- 表面反射光 --------------
        // MicrofacetDistribution distr(m_type,
        //                              m_alpha_u->eval_1(si, active),
        //                              m_alpha_v->eval_1(si, active),
        //                              m_sample_visible);
        // // Walter, et al
        // MicrofacetDistribution sample_distr(distr);
        // if (unlikely(!m_sample_visible))
        //     sample_distr.scale_alpha(1.2f - .2f * dr::sqrt(dr::abs(cos_theta_i)));

        // // Microfacet法線をサンプリング
        // Normal3f m;
        // std::tie(m, bs.pdf) =
        //     sample_distr.sample(dr::mulsign(si.wi, cos_theta_i), sample2);
        // active &= dr::neq(bs.pdf, 0.f);

        // auto [F, cos_theta_t, eta_it, eta_ti] =
        //     fresnel(dr::dot(si.wi, m), m_eta);

        // // Select the lobe to be sampled
        // UnpolarizedSpectrum weight;
        // Mask selected_r, selected_t;
        // if (likely(has_reflection && has_transmission)) {
        //     selected_r = sample1 <= F && active;
        //     weight = 1.f;
        //     /* For differentiable variants, lobe choice has to be detached to avoid bias.
        //         Sampling weights should be computed accordingly. */
        //     if constexpr (dr::is_diff_v<Float>) {
        //         if (dr::grad_enabled(F)) {
        //             weight = dr::select(selected_r, F / dr::detach(F), (1 - F) / (1.f - dr::detach(F)));
        //         }
        //     }
        //     bs.pdf *= dr::detach(dr::select(selected_r, F, 1.f - F));
        // } else {
        //     if (has_reflection || has_transmission) {
        //         selected_r = Mask(has_reflection) && active;
        //         weight = has_reflection ? F : (1.f - F);
        //     } else {
        //         return { bs, 0.f };
        //     }
        // }

        // selected_t = !selected_r && active;

        // bs.eta               = dr::select(selected_r, Float(1.f), eta_it);
        // bs.sampled_component = dr::select(selected_r, UInt32(0), UInt32(1));
        // bs.sampled_type      = dr::select(selected_r,
        //                               UInt32(+BSDFFlags::GlossyReflection),
        //                               UInt32(+BSDFFlags::GlossyTransmission));

        // Float dwh_dwo = 0.f;

        // // Reflection sampling
        // if (dr::any_or<true>(selected_r)) {
        //     // Perfect specular reflection based on the microfacet normal
        //     bs.wo[selected_r] = reflect(si.wi, m);

        //     if (m_specular_reflectance)
        //         weight[selected_r] *= m_specular_reflectance->eval(si, selected_r);

        //     // Jacobian of the half-direction mapping
        //     dwh_dwo = dr::rcp(4.f * dr::dot(bs.wo, m));
        // }

        // // Transmission sampling
        // if (dr::any_or<true>(selected_t)) {
        //     // Perfect specular transmission based on the microfacet normal
        //     bs.wo[selected_t]  = refract(si.wi, m, cos_theta_t, eta_ti);

        //     /* For transmission, radiance must be scaled to account for the solid
        //        angle compression that occurs when crossing the interface. */
        //     UnpolarizedSpectrum factor = (ctx.mode == TransportMode::Radiance) ? dr::sqr(eta_ti) : Float(1.f);

        //     if (m_specular_transmittance)
        //         factor *= m_specular_transmittance->eval(si, selected_t);

        //     weight[selected_t] *= factor;

        //     // Jacobian of the half-direction mapping
        //     dr::masked(dwh_dwo, selected_t) =
        //         (dr::sqr(bs.eta) * dr::dot(bs.wo, m)) /
        //          dr::sqr(dr::dot(si.wi, m) + bs.eta * dr::dot(bs.wo, m));
        // }

        // if (likely(m_sample_visible))
        //     weight *= distr.smith_g1(bs.wo, m);
        // else
        //     weight *= distr.G(si.wi, bs.wo, m) * dr::dot(si.wi, m) /
        //               (cos_theta_i * Frame3f::cos_theta(m));

        // bs.pdf *= dr::abs(dwh_dwo);

        // // return { bs, depolarizer<Spectrum>(weight) & active };

        
        // --------------- 再帰反射光 --------------------
        UInt32 index = (theta_i * 180.0 / dr::Pi<Float>) / 5.0f;
        index = dr::select(index >= 10, 9, index);
        index = dr::select(index < 0, 0, index);
        Float index_float = (theta_i * 180.0 / dr::Pi<Float>) / 5.0f;
        index_float = dr::select(index_float >= 10.0, 9.0, index_float);
        index_float = dr::select(index < 0.0, 0.0, index_float);
        Float index_decimal = index_float - index;

        // BSSRDF
        Point3f p = si.p; Point2f r2f = rand2(rng);
        Float sampled_p1f; Float pdfvalue; Float sampled;
        Float r1 = rand(rng); Float r2 = rand(rng);
        // sampled_p2f = m_rayshiftLUT.sample(r1);
        
        printf("cos(theta_i) = %f\n", cos_theta_i);
        printf("index = %d\n", index);

        std::tie(sampled_p1f, pdfvalue, sampled) = m_rayshiftLUT_list.sample1D(sample1, index, active);
        
        float cornersize = m_cornersize.get()->max();
        float shiftoffset = m_shiftoffset.get()->max();
        Float shift = (cornersize*sampled_p1f*1.50f)/rayshift_N; // TODO あとでちゃんとここを作る
        Float r3 = rand(rng);
        Float theta = 2*r3*dr::Pi<Float>;

        bs.p = Point3f(shift*dr::cos(theta), shift*dr::sin(theta), shiftoffset);

        // ----- 出射光のサンプリング -----
        Float point; Float sincvalue; Float sample; Float point_p1; Float sincvalue_p1; Float sample_p1;
        std::tie(point, sincvalue, sample) = m_pdfdata.sample1D(sample2.y(), index, active);
        point = dr::select(point >= 0.f, point, 0.f); // point < 0のとき0にする
        std::tie(point_p1, sincvalue_p1, sample_p1) = dr::select(index < 9, m_pdfdata.sample1D(sample, index+1, active), m_pdfdata.
        sample1D(sample, index, active));
        point_p1 = dr::select(point_p1 >= 0.f, point_p1, 0.f); // point < 0のとき0にする
        Float point_blend = point*(1-index_decimal) + point_p1*index_decimal;
        
        Float random2 = rand(rng), random3 = rand(rng);
        Float del_phi = 2.0*dr::Pi<Float> * random2; // 角度の差分にする

        Float random1 = rand(rng);
        Float offset2 = 0.0001555555555554644*dr::Pi<Float>;
        Float del_theta = dr::Pi<Float>*point_blend/(2.0*N) + offset2;

        Float del_phi_offset = (2.0f*dr::Pi<Float>)*random3; // 角度の差分にする
        Float cos_offset_phi = dr::cos(del_phi_offset), sin_offset_phi = dr::sin(del_phi_offset);
        Float cos_offset_theta = dr::cos(offset), sin_offset_theta = dr::sin(offset);
        Vector3f offsetvec = Vector3f(sin_offset_theta*cos_offset_phi, sin_offset_theta*sin_offset_phi, cos_offset_theta);        
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_theta);

        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;

        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        bs.eta = 1.f;

        Float ref = reflist.eval_pmf(index, active);
        Float ref_p1 = dr::select(index < 9, reflist.eval_pmf(index+1, active), reflist.eval_pmf(9, active));
        Float ref_blend = ref*(1-index_decimal) + ref_p1*index_decimal;
        UnpolarizedSpectrum value = Vector3f(ref_blend, ref_blend, ref_blend);
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active);
        bs.pdf = 1.f;

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

        // Float arg = (m_a->eval(si, active).x()*300/dr::cos(offset)) * dr::tan(delta_theta - offset);

        float a, reflectance;
        // if (m_isEstimation) {a = m_a.get()->max(); reflectance = m_reflectance.get()->max();}
        // else {a = alist[angle_index]; reflectance = reflist[angle_index];}
        a = m_a.get()->max(); reflectance = m_reflectance.get()->max();
        
        Float arg = (a*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
        Float sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            ScalarFloat(1.0),
            sincvalue
        );
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        // UnpolarizedSpectrum value = reflectance * sincvalue;

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

        float a, reflectance;
        // if (m_isEstimation) {a = m_a.get()->max(); reflectance = m_reflectance.get()->max();}
        // else {a = alist[angle_index]; reflectance = reflist[angle_index];}
        a = m_a.get()->max(); reflectance = m_reflectance.get()->max();

        Float arg = (a*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
        Float sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            ScalarFloat(1.0),
            sincvalue
        );
        return sincvalue;
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

            // Float arg_Float = (m_a->eval(si, active).x()*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
            
            float a, reflectance;
            // if (m_isEstimation) {a = m_a.get()->max(); reflectance = m_reflectance.get()->max();}
            // else {a = alist[angle_index]; reflectance = reflist[angle_index];}
            a = m_a.get()->max(); reflectance = m_reflectance.get()->max();

            Float arg_Float = (a*300/dr::cos(offset)) * dr::tan(delta_theta - offset);
            Float sincvalue_pdf = dr::pow(dr::sin(arg_Float), 2) / dr::pow(arg_Float, 2);
            sincvalue_pdf = dr::select(
                dr::isnan(sincvalue_pdf),
                ScalarFloat(1.0),
                sincvalue_pdf
            );
            UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue_pdf;
            // UnpolarizedSpectrum value = m_reflectance->eval(si, active);
            // UnpolarizedSpectrum value = reflectance;
            Float pdf = 1.f;

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
        }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "sinc_shift4[" << std::endl
            << " reflectance = " << string::indent(m_reflectance) << std::endl
            << "a = " << string::indent(m_a) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

// フィールド
private:
    // 表面反射光関係
    ref<Texture> m_specular_reflectance;
    ref<Texture> m_specular_transmittance;
    MicrofacetType m_type;
    ref<Texture> m_alpha_u, m_alpha_v;
    Float m_eta, m_inv_eta;
    bool m_sample_visible;

    // 再帰反射光関係
    uint32_t M; // LUTのphiサイズ
    uint32_t N; // LUTのthetaサイズ
    uint32_t rayshift_N; // rayshift_LUTのphiサイズ
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    float a, reflectance;
    int angle, angle_index;
    ScalarFloat *m_data;
    DiscreteDistribution2D<Float, 2> m_pdfdata;
    DiscreteDistribution2D<Float, 2> m_rayshiftLUT_list;
    ScalarFloat offset;
    ref<Texture> m_angle;
    ref<Texture> m_cornersize;
    ref<Texture> m_shiftoffset;
    bool m_isEstimation;
    float alist[10] = {8.97515869140625, 10.97819995880127, 10.85496711730957, 8.351812362670898, 8.03887939453125, 4.881270408630371, 3.225633144378662, 2.2548987865448, 1.620800495147705};
    DiscreteDistribution<Float> reflist;
};

MI_IMPLEMENT_CLASS_VARIANT(sinc_shift4, BSDF)
MI_EXPORT_PLUGIN(sinc_shift4, "sinc_shift4")
NAMESPACE_END(mitsuba)