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

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <random>
#include <iterator>
#include <string>
#include "nlohmann/json.hpp"
using namespace std;
using json = nlohmann::json;


NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class sinc_shift3 final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    using Index = dr::uint32_array_t<Float>;
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    sinc_shift3(const Properties &props) : Base(props)
    {
        m_rrrw = m_rrrh = 175.0;
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        m_divergence = props.texture<Texture>("divergence", 1.f);
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        m_a = props.texture<Texture>("a", 1.f);
        m_angle = props.texture<Texture>("angle", 0.f);
        m_cornersize = props.texture<Texture>("cornersize", 0.f);
        m_shiftoffset = props.texture<Texture>("shiftoffset", 0.f);

        M = 1; N = 30000;

        // offset = dr::Pi<Float>/15555.0f;
        // offset = 0.0;
        // offset = 0.0001000*dr::Pi<Float>;
        // offset = 0.00012500*dr::Pi<Float>;        
        offset = 0.0001555555555554644*dr::Pi<Float>;
        // offset = 18.0/(3600.0*180.0)*dr::Pi<Float>;
        // offset = 0.00010*dr::Pi<Float>;
        // offset = 0.00019055555555558268*dr::Pi<ScalarFloat>;
        // offset = 0.0002*dr::Pi<Float>;
        // offset = 0.00025*dr::Pi<Float>;
        // offset = 0.0003*dr::Pi<Float>;
        // offset = 0.0004*dr::Pi<Float>;
        // offset = 0.0005*dr::Pi<Float>;
        // offset = 0.0006*dr::Pi<Float>;
        // offset = 0.0008*dr::Pi<Float>;
        // offset = dr::Pi<Float>/15555.0f;
        // offset = 0.01*dr::Pi<Float>;
        // offset = 0.f;

        make_rayshift_array();

        isTraversed = false;
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
        ScalarFloat data[N];
        int i = 0, idx = 0;

        for (int i = 0; i < N; ++i)
        {
            ScalarFloat arg = a * 300*dr::tan((dr::Pi<Float>*i/(2.0f*N) - offset));
            // ScalarFloat arg = a * dr::tan(dr::Pi<Float>*i/(2.5f*N) - offset);
            // ScalarFloat arg = a * 300*dr::tan((dr::Pi<Float>*i/(2.0f*N)));
            // ScalarFloat arg = a *dr::tan((dr::Pi<Float>*i/(2.0f*N)));
            // ScalarFloat arg = a * (-dr::Pi<Float>/2.0f + dr::Pi<Float>*i/N);
            // ScalarFloat arg = a * (dr::Pi<Float>*i/(10.0*N));
            // ScalarFloat arg = a * (20.0f*i/N - 10.0f);
            ScalarFloat sincvalue = dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo);
            // sincvalue = dr::select(
            //     dr::isfinite(sincvalue),
            //     dr::abs(sincvalue),
            //     ScalarFloat(1.0)
            // );
            if (i >= 28000) sincvalue = 0.f;
            // if (i >= N-10) printf("%f\n", sincvalue);
            // printf("%d = %f ", i, sincvalue);
            sincvalue = dr::select(
                dr::isinf(sincvalue),
                ScalarFloat(1.0),
                dr::abs(sincvalue)
            );
            data[i] = sincvalue;
        }
        // 正規化（ごり押し）
        float sum = 0.0f;
        for (int i = 0; i < N; ++i)
        {
            sum += data[i];
        }
        for (int i = 0; i < N; ++i)
        {
            data[i] /= sum;
        }

        sum_rawPDF = sum;
        m_data = &data[0]; // LUTの先頭ポインタ
        struct DiscreteDistribution<Float> dd(m_data, N);
        m_pdfdata = dd;
    }

    void make_rayshift_array() {
        std::string filename = "/home/sugawara.ryo/makeLUT/2025.03.17diff/LUTs/rayshift_LUT.json";
        ifstream ifs(filename.c_str());
        if (ifs.good())
        {
            json m_json;
            ifs >> m_json;

            // for (const auto& json : m_json) {
            //     // vector<vector<float>> json;
            //     for (const auto& each_vector : json)
            //     {
            //         vector<float> d = each_vector.get<vector<float>>();
            //         m_rayshift_list.push_back(d);
            //     }
            // }

            for (const auto& json : m_json) {
                vector<float> d = json.get<vector<float>>();
                m_rayshift_list.push_back(d);
            }
        }
        else
        {
            cout << "ファイルの読み込みに失敗しました" << endl;
        }

        int num_i = m_rayshift_list.size();
        int num_j = m_rayshift_list[0].size();
        // ScalarFloat data[num_i*num_j];
        ScalarFloat data[num_j];
        int i = 0, idx = 0;
        // for (int i = 0; i < num_i; ++i) {
        //     for (int j = 0; j < num_j; ++j) {
        //         data[idx++] = m_rayshift_list[i][j];
        //     }
        // }
        for (int j = 0; j < num_j; ++j) {
            data[idx++] = m_rayshift_list[9][j];
        }

        // for (int i = 0; i < num_i*num_j; ++i) printf("data[idx] = %f\n", data[i]);

        // ScalarVector2u size(num_i, num_j);
        ScalarVector2u size(num_j, num_i);
        struct DiscreteDistribution<Float> dd(&data[0], num_j);
        // struct DiscreteDistribution2D<Float, 2> dd(&data[0], size);
        m_rayshiftLUT = dd;
        rayshift_N = num_j;
        printf("rayshift_N = %d\n", num_j);
    }

    void traverse(TraversalCallback *callback) override
    {
        callback->put_object("divergence", m_divergence.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("angle", m_angle.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("cornersize", m_cornersize.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("shiftoffset", m_shiftoffset.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);

        isTraversed = true;

        float a = m_a.get()->max();
        int angle = static_cast<int>(m_angle.get()->max());        
        
        make_sincarray1(m_a.get()->max(), 2, angle);
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

        // 光線がぶつかった点からコーナーキューブのどの点に入射したかを判定
        // Point3f p = si.p; Float corner_i = m_rrrw / (cornersize*3.0/2.0);

        // BSSRDF
        Point3f p = si.p; Point2f r2f = rand2(rng);
        Float sampled_p2f; Float pdfvalue; Point2f sampled;
        Float r1 = rand(rng); Float r2 = rand(rng);
        sampled_p2f = m_rayshiftLUT.sample(r1);

        float cornersize = m_cornersize.get()->max();
        float shiftoffset = m_shiftoffset.get()->max();
        Float shift = (cornersize*sampled_p2f*1.50f)/rayshift_N; // TODO あとでちゃんとここを作る
        Float r3 = rand(rng);
        Float theta = 2*r3*dr::Pi<Float>;

        bs.p = Point3f(shift*dr::cos(theta), shift*dr::sin(theta), shiftoffset);

        // ----- 出射光のサンプリング -----
        Float point;
        point = m_pdfdata.sample(sample2.y());
        Float random2 = rand(rng), random3 = rand(rng);
        Float del_phi = 2.0*dr::Pi<Float> * random2; // 角度の差分にする
        // Float sign = dr::select(random3 >= 0.5, 1.0f, -1.0f);
        // Float del_phi_offset = theta; // 角度の差分にする
        // Float del_phi_offset = (2.0f*dr::Pi<Float>)*(( ((100*random3) % 6) + 1) / 6 ); // 角度の差分にする
        // Float a = m_a->eval(si, active).x();
        Float random1 = rand(rng);
        // thetaサンプリングバージョン
        // Float del_theta = point*dr::Pi<Float> / N - dr::Pi<Float>/2.0f;
        // Float del_theta = dr::Pi<Float>*point/(10.0*N) + offset;
        // Float del_theta = dr::Pi<Float>*point/(2.5*N) + offset;
        // printf("point = %d\n", point);
        Float del_theta = dr::Pi<Float>*point/(2.0*N);
        // printf("del_theta = %f", del_theta);
        
        // xサンプリングバージョン
        // Float X = 20.0f*point/N - 10.0f;
        // Float x = 30000.0f*point/N - 15000.0f;
        // Float x_offset = shift - cornersize;
        // Float del_theta = dr::atan(X*dr::cos(offset)/300.0f) + offset;
        // Float del_theta = dr::atan(X/300.0f);
        // Float del_theta = dr::atan((2*x - (x_offset)*dr::cos(offset)) / (2*300*dr::cos(offset) - x_offset*sin(offset)));
        // Float del_theta = dr::atan((2*x - (x_offset)) / (2*300));
        // Float del_theta = dr::atan((x*dr::cos(offset))/300.0) + offset;
        // Float del_theta = point*dr::Pi<Float>/N - dr::Pi<Float>/2.0f;
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N);

        // printf("x = %f, del_theta = %f\n", x, del_theta);

        Float del_phi_offset = (2.0f*dr::Pi<Float>)*random3; // 角度の差分にする
        Float cos_offset_phi = dr::cos(del_phi_offset), sin_offset_phi = dr::sin(del_phi_offset);
        Float cos_offset_theta = dr::cos(offset), sin_offset_theta = dr::sin(offset);
        Vector3f offsetvec = Vector3f(sin_offset_theta*cos_offset_phi, sin_offset_theta*sin_offset_phi, cos_offset_theta);        
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_theta);

        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        // Vector3f wo = rotate(delvec, norm, offsetvec);
        // wo = rotate(wo, norm, si.wi);
        Vector3f wo = rotate(delvec, norm, si.wi);
        
        bs.wo = wo;

        // bs.sampled_component = 0;
        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        bs.eta = 1.f;
        // bs.pdf = 1.f;

        // UnpolarizedSpectrum arg = m_a->eval(si, active)*(20.0f*point/N - 10.0f);
        UnpolarizedSpectrum arg = m_a->eval(si, active)*300*dr::tan((del_theta - offset));
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*dr::tan((del_theta - offset));
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // UnpolarizedSpectrum sincvalue = dr::select(arg != 0.f, dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2), 1.f);
        Float arg_pdf = m_a->eval(si, active).x()*300*dr::tan((del_theta - offset));
        // Float arg_pdf = m_a->eval(si, active).x()*dr::tan((del_theta - offset));
        Float sincvalue_pdf = dr::pow(dr::sin(arg_pdf), 2) / dr::pow(arg_pdf, 2);
        // Float sincvalue_pdf = dr::select(arg_pdf != 0.f, dr::pow(dr::sin(arg_pdf), 2) / dr::pow(arg_pdf, 2), 1.f);
        // printf("pdf = %f\n", bs.pdf);
        // UnpolarizedSpectrum sincvalue = 1.0;
        // Mask iszero = (arg == 0.0);
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            // iszero,
            sincvalue,
            ScalarFloat(1.0)
        );
        sincvalue_pdf = dr::select(
            // dr::isfinite(sincvalue_pdf),
            dr::isinf(sincvalue_pdf),
            // iszero,
            // sincvalue_pdf,
            ScalarFloat(1.0),
            sincvalue_pdf
            // ScalarFloat(1.0)
        );
        // printf("arg = %f, sincvalue = %f\n", arg_pdf, sincvalue_pdf);

        // sincvalue = 1.0;

        // sincvalue_pdf = 1.f;
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue_pdf;
        UnpolarizedSpectrum value = m_reflectance->eval(si, active);
        // bs.pdf = sincvalue_pdf;
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
        Float delta_theta = acos(cos_delta_theta) - offset;

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        Float a = m_a->eval(si, active).x();
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*(300*dr::tan(delta_theta - offset));
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*(300*dr::tan(delta_theta));
        Float arg = m_a->eval(si, active).x()*(300*dr::tan(delta_theta - offset));
        // Float arg = m_a->eval(si, active).x()*(300*dr::tan(delta_theta));
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*(dr::tan(delta_theta));
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta);
        // UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        Float sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask iszero = (arg == 0.0);
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            // iszero,
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue = 1.0;
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue * 300.0 / cos_delta_theta;
        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;

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
        Float delta_theta = acos(cos_delta_theta) - offset;

        Float a = m_a->eval(si, active).x();
        Float arg = a*(300*dr::tan(delta_theta - offset));
        // Float arg = a*(300*dr::tan(delta_theta));
        // Float arg = a*(dr::tan(delta_theta));
        // Float arg = a*(delta_theta);
        Float sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask iszero = (arg == 0.0);
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            // iszero,
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue = 1.0;

        // return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue.x(), 0.f);
        Float cos_theta_del = cos_theta_i*cos_theta_o + sin_theta_i*sin_theta_o;
        // return sincvalue * 300.0/(cos_delta_theta*cos_delta_theta);
        return sincvalue;
        // return 1.f;
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
            Float delta_theta = acos(cos_delta_theta) - offset;

            active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

            Float a = m_a->eval(si, active).x();
            // UnpolarizedSpectrum arg = m_a->eval(si, active)*(300*dr::tan(delta_theta - offset));
            UnpolarizedSpectrum arg = m_a->eval(si, active)*(300*dr::tan(delta_theta));
            // UnpolarizedSpectrum arg = m_a->eval(si, active)*(dr::tan(delta_theta));
            // UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta);
            UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
            Float arg_Float = m_a->eval(si, active).x()*(300*dr::tan(delta_theta - offset));
            // Float arg_Float = m_a->eval(si, active).x()*(dr::tan(delta_theta));
            // Float arg_Float = m_a->eval(si, active).x()*(delta_theta);
            Float sincvalue_pdf = dr::pow(dr::sin(arg_Float), 2) / dr::pow(arg_Float, 2);
            // Float sincvalue_pdf = 1.0f;
            // Mask iszero = (arg == 0.0);
            // printf("arg = %f, sincvalue = %f\n", dr::tan(delta_theta - offset), sincvalue);
            sincvalue = dr::select(
                dr::isfinite(sincvalue),
                // iszero,
                sincvalue,
                ScalarFloat(1.0)
            );
            sincvalue_pdf = dr::select(
                dr::isfinite(sincvalue_pdf),
                // iszero,
                sincvalue_pdf,
                ScalarFloat(1.0)
            );
            // sincvalue_pdf = 1.f;
            // sincvalue = 1.0;
            // UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue * 300.0 / cos_delta_theta;
            // UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue_pdf;
            UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue_pdf;
            // Float pdf = dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue.x(), 0.f);
            // Float pdf = sincvalue_pdf * 300.0 / (cos_delta_theta*cos_delta_theta);
            Float pdf = sincvalue_pdf;
            // Float pdf = 1.f;

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
        }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "sinc_shift3[" << std::endl
            << " divergence = " << string::indent(m_divergence) << std::endl
            << " reflectance = " << string::indent(m_reflectance) << std::endl
            << "a = " << string::indent(m_a) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

// フィールド
private:
    uint32_t M; // LUTのphiサイズ
    uint32_t N; // LUTのthetaサイズ
    uint32_t rayshift_N; // rayshift_LUTのphiサイズ
    float sum_rawPDF;
    ref<Texture> m_divergence;
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    ScalarFloat *m_data;
    DiscreteDistribution<Float> m_pdfdata;
    vector<vector<float>> m_rayshift_list;
    // DiscreteDistribution2D<Float, 2> m_rayshiftLUT;
    DiscreteDistribution<Float> m_rayshiftLUT;
    std::string LUTFilename;
    std::string LUTDir;
    bool isTraversed;
    ScalarFloat offset;
    Float LUTWidth;
    ref<Texture> m_angle;
    ref<Texture> m_cornersize;
    ref<Texture> m_shiftoffset;
    Float m_rrrw;
    Float m_rrrh;
};

MI_IMPLEMENT_CLASS_VARIANT(sinc_shift3, BSDF)
MI_EXPORT_PLUGIN(sinc_shift3, "sinc_shift3")
NAMESPACE_END(mitsuba)