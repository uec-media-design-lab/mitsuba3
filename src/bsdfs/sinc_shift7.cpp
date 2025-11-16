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
#include <cmath>
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
    MI_IMPORT_TYPES(Texture, MicrofacetDistribution)

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
        R = 600.f;

        m_isEstimation = 1;
        // theta_max = dr::Pi<Float>/2.f;

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
        ScalarFloat sums_local[10];

        for (int i = 0; i < 10; ++i)
        {
            if (!m_isEstimation)
                a = alist[i];
            
            sums_local[i] = 0.0f;
            for (int j = 0; j < N; ++j) {
                ScalarFloat theta = dr::Pi<Float>*j/(2.0f*N);
                // ScalarFloat arg = a * dr::tan(theta);
                ScalarFloat arg = a * R * dr::tan(theta);
                // ScalarFloat sincvalue = dr::pow(dr::sin(a*dr::tan(theta)), 2) / (dr::sin(theta)*dr::cos(theta));
                // ScalarFloat sincvalue = dr::pow(dr::sin(a*dr::tan(theta)), 2) / (a*a*dr::sin(theta)*dr::sin(theta)*dr::cos(theta));
                // ScalarFloat sincvalue = dr::pow(dr::sin(a*dr::tan(theta)), 2) / (dr::pow(a*dr::tan(theta), 2)*dr::pow(dr::cos(theta), 2));

                // 角度ベースサンプリング
                // ScalarFloat sincvalue = dr::pow(dr::sin(a*dr::tan(theta)), 2) / (a*a*dr::sin(theta)*dr::cos(theta));
                ScalarFloat sincvalue = dr::pow(dr::sin(arg), 2) / (a*a*dr::sin(theta)*dr::cos(theta));
                // ScalarFloat sincvalue = dr::pow(dr::sin(arg), 2) / (a*a*dr::sin(theta)*dr::sin(theta)*dr::cos(theta));

                // sincvalue = dr::select(
                //     dr::isfinite(sincvalue),
                //     dr::abs(sincvalue),
                //     ScalarFloat(1.0)
                // );
                if (j == 0) sincvalue = 0.0f;
                // if (theta > dr::atan(dr::Pi<Float>/3.f)) sincvalue = 0.f;
                // if (theta > dr::atan(10.f*dr::Pi<Float>/(a*R))) sincvalue = 0.f;
                if (theta > dr::atan(2.4f*dr::Pi<Float>/(a*R))) sincvalue = 0.f;
                // if (theta > 1.8185301*dr::Pi<Float>/180.f) sincvalue = 0.f;
                // if (j >= 28000) sincvalue = 0.f;
                // if (j == N-1) sincvalue = 0.0f;
                data[idx++] = sincvalue;
                sums_local[i] += sincvalue;
            }
        }

        // 正規化（ごり押し）
        idx = 0;
        for (int i = 0; i < 10; ++i) {
            // printf("sums_local[%d] = %f\n", i, sums_local[i]);
            for (int j = 0; j < N; ++j)
            {   
                // if (idx < N) printf("data[%d] = %f\n", idx, data[idx]);
                data[idx++] /= sums_local[i];
            }
        }

        DiscreteDistribution2D<Float, 2> dd(&data[0], size);
        m_pdfdata = dd;

        float reflist_float[10] = {0.251913845539093, 0.24916136264801025, 0.24612076580524445, 0.2393658608198166, 0.23089219629764557, 0.22330009937286377, 0.20480982959270477, 0.18945062160491943, 0.16936106979846954, 0.2019372433423996};

        ScalarFloat data2[10];
        for (int i = 0; i < 10; ++i)
            data2[i] = reflist_float[i];
        DiscreteDistribution<Float> rl(&data2[0], 10);
        DiscreteDistribution<Float> sums_dist(&sums_local[0], 10);
        sums = sums_dist;
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
                data[idx++] = m_rayshift_list[i][j];
            }
        }
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

        bool has_reflection    = ctx.is_enabled(BSDFFlags::GlossyReflection, 0),
        has_transmission  = ctx.is_enabled(BSDFFlags::GlossyTransmission, 1);

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};

        auto rng = setRandomGenerator(sample1*10000000);

        // ----------- 表面反射光 --------------
        MicrofacetDistribution distr(m_type,
                                     m_alpha_u->eval_1(si, true),
                                     m_alpha_v->eval_1(si, true),
                                     m_sample_visible);
        MicrofacetDistribution sample_distr(distr);
        
        // Walter, et al
        if (unlikely(!m_sample_visible))
            sample_distr.scale_alpha(1.2f - .2f * dr::sqrt(dr::abs(cos_theta_i)));

        // Microfacet法線をサンプリング
        Normal3f m; Float Di = 1.f, Di_in = 1.f;
        // printf("Di-0 = %f\n", Di);
        // Di_inはWalterの論文の中のp_m
        std::tie(m, Di_in) = sample_distr.sample(dr::mulsign(si.wi, cos_theta_i), sample2);
        active &= dr::neq(Di_in, 0.f); Di *= Di_in;
        // printf("Di-1 = %f\n", Di);
        auto [F, cos_theta_t, eta_it, eta_ti] = fresnel(dr::dot(si.wi, m), m_eta); // フレネル項

        // Select the lobe to be sampled
        UnpolarizedSpectrum weight;
        Mask selected_r, selected_t;
        if (likely(has_reflection && has_transmission)) {
            selected_r = sample1 <= F;
            weight = 1.f;
            if constexpr (dr::is_diff_v<Float>) {
                if (dr::grad_enabled(F)) {
                    weight *= dr::select(selected_r, F / dr::detach(F), (1 - F) / (1.f - dr::detach(F)));
                }
            }
            Di *= dr::detach(dr::select(selected_r, F, 1.f - F));
        } else {
            if (has_reflection || has_transmission) {
                selected_r = Mask(has_reflection) && active;
                weight = has_reflection ? F : (1.f - F);
            } else {
                return { bs, 0.f };
            }
        }
        // printf("Di-2 = %f\n", Di);

        // 表面を透過する場合
        selected_t = !selected_r;
        bs.eta               = dr::select(selected_r, Float(1.f), eta_it);
        bs.sampled_component = dr::select(selected_r, UInt32(0), UInt32(1));
        bs.sampled_type      = dr::select(selected_r,
                                      UInt32(+BSDFFlags::GlossyReflection),
                                      UInt32(+BSDFFlags::GlossyTransmission));
        
        Float dwh_dwo = 0.f; // Jacobian

        Vector3f norm = si.to_local(si.n);
        // 表面反射する場合
        if (dr::any_or<true>(selected_r)) {
            bs.wo[selected_r] = reflect(si.wi, m);
            if (m_specular_reflectance)
                weight[selected_r] *= m_specular_reflectance->eval(si, selected_r);

            // Jacobian of the half-direction mapping
            dwh_dwo = dr::rcp(4.f * dr::dot(bs.wo, m));
            Di *= dr::abs(dwh_dwo);
            // printf("Di-3-ref = %f\n", Di);
            bs.p[selected_r] = Point3f(0.f, 0.f, 0.f);
            if (likely(m_sample_visible))
                weight[selected_r] *= distr.smith_g1(bs.wo, m);
            else
                weight[selected_r] *= distr.G(si.wi, bs.wo, m) * dr::dot(si.wi, m) /
                        (cos_theta_i * Frame3f::cos_theta(m));
        }

        // 透過する場合
        if (dr::any_or<true>(selected_t)) {
            // 最初の屈折
            Vector3f v1 = refract(si.wi, m, cos_theta_t, eta_ti);
            UnpolarizedSpectrum factor = (ctx.mode == TransportMode::Radiance) ? dr::sqr(eta_ti) : Float(1.f);
            if (m_specular_transmittance)
                factor *= m_specular_transmittance->eval(si, selected_t);
            weight[selected_t] *= factor;
            dr::masked(dwh_dwo, selected_t) =
                (dr::sqr(bs.eta) * dr::dot(v1, m)) /
                 dr::sqr(dr::dot(si.wi, m) + bs.eta * dr::dot(v1, m));
            Di *= dr::select(selected_t, dr::abs(dwh_dwo), 1.f);
            // printf("Di-3-trans = %f\n", Di);
            if (likely(m_sample_visible))
                weight[selected_t] *= distr.smith_g1(bs.wo, m);
            else
                weight[selected_t] *= distr.G(si.wi, v1, m) * dr::dot(si.wi, m) /
                        (cos_theta_i * Frame3f::cos_theta(m));
            // printf("weight_1 = %.15f\n", weight);
            
            // 再帰反射光
            // 入射角は0 - 45度のうちどこか
            UInt32 index = dr::floor((theta_i * 180.0 / dr::Pi<Float>) / 5.0f);
            index = dr::select(index >= 10, 9, index); index = dr::select(index < 0, 0, index);
            Float index_float = (theta_i * 180.0 / dr::Pi<Float>) / 5.0f;
            index_float = dr::select(index_float >= 10.0, 9.0, index_float);
            index_float = dr::select(index < 0.0, 0.0, index_float);
            Float index_decimal = index_float - index;

            // BSSRDF
            Point3f p = si.p; Point2f r2f = rand2(rng);
            Float sampled_p1f; Float pdfvalue; Float sampled;
            Float r1 = rand(rng); Float r2 = rand(rng);
            std::tie(sampled_p1f, pdfvalue, sampled) = m_rayshiftLUT_list.sample1D(sample1, index, selected_t);
            float cornersize = m_cornersize.get()->max();
            float shiftoffset = m_shiftoffset.get()->max();
            Float shift = (cornersize*sampled_p1f*1.50f)/rayshift_N; // TODO あとでちゃんとここを作る
            Float r3 = rand(rng);
            Float theta = 2*r3*dr::Pi<Float>;
            bs.p = Point3f(shift*dr::cos(theta), shift*dr::sin(theta), shiftoffset);

            // 出射方向サンプリング
            Float point; Float sincvalue; Float sample; Float point_p1; Float sincvalue_p1; Float sample_p1;
            std::tie(point, sincvalue, sample) = dr::select(index <= 9, m_pdfdata.sample1D(sample2.y(), index, selected_t), m_pdfdata.sample1D(sample2.y(), 9, selected_t));
            point = dr::select(point >= 0.f, point, 0.f); // point < 0のとき0にする

            Float random2 = rand(rng), random3 = rand(rng);
            Float del_phi = 2.0*dr::Pi<Float> * random2; // 角度の差分にする

            Float random1 = rand(rng);
            Float offset2 = 0.f;
            Float del_theta = dr::Pi<Float>*point/(2.0f*N) + offset2;
            // Float del_theta = 0.f;
            Float a_Float = m_a->eval(si, selected_t).x();
            // Float arg = a_Float*dr::tan(del_theta);
            Float arg = a_Float*R*dr::tan(del_theta);
            sincvalue = dr::pow(dr::sin(arg), 2)/( a_Float*a_Float* dr::cos(del_theta)*dr::pow(dr::sin(del_theta), 2));
            // sincvalue = dr::pow(dr::sin(arg), 2)/(dr::pow(arg, 2)*dr::pow(dr::cos(del_theta), 2));
            sincvalue /= sums.eval_pmf(index, selected_t);
            Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
            Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
            Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_theta);
            Vector3f v2 = rotate(delvec, norm, v1); // 再帰反射方向

            if (m_isEstimation) {
                // weight[selected_t] = m_reflectance->eval(si, selected_t)*sincvalue;
                // weight[selected_t] = m_reflectance->eval(si, selected_t)*Frame3f::cos_theta(v1);
                weight[selected_t] = m_reflectance->eval(si, selected_t);
            } else {
                Float ref = dr::select(index <= 9, reflist.eval_pmf(index, selected_t), reflist.eval_pmf(9, selected_t));
                Float ref_p1 = dr::select(index < 9, reflist.eval_pmf(index+1, selected_t), reflist.eval_pmf(9, selected_t));
                Float ref_blend = ref*(1-index_decimal) + ref_p1*index_decimal;
                // sincvalue *= sums.eval_pmf(index, selected_t);
                // weight[selected_t] = Vector3f(ref_blend, ref_blend, ref_blend)*sincvalue*Frame3f::cos_theta(v1);
                weight[selected_t] = Vector3f(ref, ref, ref)*Frame3f::cos_theta(v1);
            }
            Di *= sincvalue;
            // Di *= 1.f;

            // 最後の透過
            auto [F, cos_theta_to, eta_ot, eta_to] = fresnel(dr::dot(-v2, m), 1.f/m_eta); // フレネル項
            UnpolarizedSpectrum factor2 = (ctx.mode == TransportMode::Radiance) ? dr::sqr(eta_ot) : Float(1.f);
            if (m_specular_transmittance)
                factor2 *= m_specular_transmittance->eval(si, selected_t);
            weight[selected_t] *= factor2;
            Vector3f wo = refract(v2, -m, cos_theta_to, eta_to);
            bs.wo[selected_t] = wo;
            dr::masked(dwh_dwo, selected_t) = 
                (dr::sqr(eta_ot) * dr::dot(wo, m)) /
                 dr::sqr(dr::dot(v2, m) + eta_ot * dr::dot(wo, m));
            Di *= dr::detach(1.f - F)*Di_in;
            // printf("Di-4 = %f\n", Di);
            Di *= dr::select(selected_t, dr::abs(dwh_dwo), 1.f);
            // printf("Di-5 = %f\n", Di);
            if (likely(m_sample_visible))
                weight[selected_t] *= distr.smith_g1(bs.wo, m);
            else
                weight[selected_t] *= distr.G(v2, bs.wo, m) * dr::dot(v2, m) /
                        (cos_theta_t * Frame3f::cos_theta(m));
            // printf("weight_2 = %f %f %f\n", weight.x(), weight.y(), weight.z());
            bs.sampled_type =+ BSDFFlags::DiffuseReflection;
            bs.eta = 1.f;
        }

        bs.pdf = Di;
        // printf("bs.pdf = %f\n", bs.pdf);

        // return {bs, depolarizer<Spectrum>(weight) & (active && bs.pdf > 0.f)};
        return {bs, depolarizer<Spectrum>(weight) & (bs.pdf > 0.f)};
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

        float a, reflectance;
        a = m_a.get()->max(); reflectance = m_reflectance.get()->max();
        
        Float arg = (a * dr::tan(delta_theta));
        Float sincvalue = dr::pow(dr::sin(arg), 2) / (a*a*dr::sin(delta_theta)*dr::cos(delta_theta));
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            ScalarFloat(1.0),
            sincvalue
        );
        UInt32 index = dr::floor((delta_theta * 180.0 / dr::Pi<Float>) / 5.0f);
        index = dr::select(index >= 10, 9, index); index = dr::select(index < 0, 0, index);
        UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue / sums.eval_pmf(index, active);
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue;

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
        a = m_a.get()->max(); reflectance = m_reflectance.get()->max();

        Float arg = (a * dr::tan(delta_theta));
        Float sincvalue = dr::pow(dr::sin(arg), 2) / (a*a*dr::sin(delta_theta)*dr::cos(delta_theta));
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            ScalarFloat(1.0),
            sincvalue
        );
        UInt32 index = dr::floor((delta_theta * 180.0 / dr::Pi<Float>) / 5.0f);
        index = dr::select(index >= 10, 9, index); index = dr::select(index < 0, 0, index);
        
        return sincvalue / sums.eval_pmf(index, active);
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
            
            float a, reflectance;
            a = m_a.get()->max(); reflectance = m_reflectance.get()->max();

            Float arg_Float = a * dr::tan(delta_theta);
            Float sincvalue_pdf = dr::pow(dr::sin(arg_Float), 2) / (a*a*dr::sin(delta_theta)*dr::cos(delta_theta));
            sincvalue_pdf = dr::select(
                dr::isnan(sincvalue_pdf),
                ScalarFloat(1.0),
                sincvalue_pdf
            );
            UnpolarizedSpectrum value = m_reflectance->eval(si, true);

            UInt32 index = dr::floor((delta_theta * 180.0 / dr::Pi<Float>) / 5.0f);
            index = dr::select(index >= 10, 9, index); index = dr::select(index < 0, 0, index);
            Float pdf = sincvalue_pdf / sums.eval_pmf(index, active);

            return {depolarizer<Spectrum>(value), pdf};
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
    float R; // 距離
    float l = 0.00000066f; // 波長
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
    DiscreteDistribution<Float> sums;
    float alist[10] = {1.2247763872146606, 1.174414038658142, 1.0820084810256958, 0.960394024848938, 0.8475515842437744, 0.6400105953216553, 0.5727102160453796, 0.41668373346328735, 0.396972626447677612, 0.3177734911441802979};
    DiscreteDistribution<Float> reflist;
    float theta_max = 3.1415926535897932384f/2.f;
};

MI_IMPLEMENT_CLASS_VARIANT(sinc_shift4, BSDF)
MI_EXPORT_PLUGIN(sinc_shift4, "sinc_shift4")
NAMESPACE_END(mitsuba)