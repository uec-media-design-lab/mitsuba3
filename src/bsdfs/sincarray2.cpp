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
// #include "nlohmann/json.hpp"
using namespace std;
// using json = nlohmann::json;


NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class sincarray2 final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    using Index = dr::uint32_array_t<Float>;
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    sincarray2(const Properties &props) : Base(props)
    {
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        m_divergence = props.texture<Texture>("divergence", 1.f);
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        m_a = props.texture<Texture>("a", 1.f);
        m_angle = props.texture<Texture>("angle", 0.f);
        m_cornersize = props.texture<Texture>("cornersize", 0.f);

        M = 1; N = 30000;

        // offset = dr::Pi<Float>/15555.0f;
        offset = 0.f;

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

    vector<vector<float>> make_sincarray2_before(vector<vector<float>>& vecs, int M, int N, float a, int expo) const
    {
        for (int i = 0; i < N; i++) { 
            vector<float> parts = vector<float>(M); // 0で初期化
            for (int j = 0; j < M; j++)
            {
                float arg = a * static_cast<float>(std::sqrt( std::pow(j - M/2, 2) + std::pow(i - N/2, 2) ));
                parts[j] = arg == 0 ? 1.0f : std::pow(std::sin(arg), expo) / std::pow(arg, expo); // argが0の時はゼロ除算を防ぐため
                cout << "parts["<< j << "] = " << parts[j] << endl;
            }
            vecs.push_back(parts);
        }
        //正規化
        vector<vector<float>> v;
        // v = normalize(vecs);
        for (int i = 0; i < N; i++) { 
            for (int j = 0; j < M; j++)
            {
                cout << "parts[" << i << "][" << j << "] = " << v[i][j] << endl;
            }
        }
        return v;
    }

    void make_sincarray2(float a, int expo)
    {
        ScalarVector2u size(M, N); ScalarFloat data[M*N];
        int i = 0, idx = 0;

        // std::cout << "a = " << m_a.get()->max() << std::endl;
        for (int i = 0; i < N; ++i)
        {
            idx = M*i;
            for (int j = 0; j < M; ++j)
            {
                // ScalarFloat arg = m_a * (dr::sqrt( dr::pow(ScalarFloat(j - M / 2.0f), 2) + dr::pow(ScalarFloat(i - N / 2.0f), 2) ));
                // ScalarFloat arg = a * (dr::sqrt( dr::pow(j - M / 2.0f, 2) + dr::pow(i - N / 2.0f, 2) ) );
                ScalarFloat arg = a * (dr::sqrt( dr::pow(dr::TwoPi<Float>*j / M, 2) + dr::pow(dr::Pi<Float>*i / (2.0f*N), 2) ) );
                ScalarFloat sincvalue = dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo);
                // Mask mask = arg == 0;
                // argが0の時はゼロ除算を防ぐため
                sincvalue = dr::select(
                    dr::isfinite(sincvalue),
                    sincvalue,
                    ScalarFloat(1.0)
                );

                data[idx] = sincvalue;
                ++idx;
            }
        }

        // 正規化（ごり押し）
        float sum = 0.0f;
        for (int i = 0; i < N; ++i)
        {
            idx = M*i;
            for (int j = 0; j < M; ++j)
            {
                sum += data[idx++];
            }
        }
        // std::cout << "[" << std::endl;
        for (int i = 0; i < N; ++i)
        {
            // std::cout << "    [" << std::endl;
            idx = M*i;
            for (int j = 0; j < M; ++j)
            {
                data[idx] /= sum;

                // if (j != M-1) std::cout << "        " << data[idx] << "," << std::endl;
                // else std::cout << "        " << data[idx] << std::endl;
                ++idx;
            }
            // if (i != N-1) std::cout << "    ]," << std::endl;
            // else std::cout << "    ]" << std::endl;
        }
        // std::cout << "]" << std::endl;

        sum_rawPDF = sum;
        m_data = &data[0]; // LUTの先頭ポインタ
        m_size = size; // LUTの縦、横サイズ
        m_pdfdata2 = DiscreteDistribution2D<Float, 2>(m_data, m_size);
    }

    void make_sincarray1(float a, int expo, int angle)
    {
        ScalarFloat data[N];
        int i = 0, idx = 0;

        for (int i = 0; i < N; ++i)
        {
            // ScalarFloat arg = a * (-dr::Pi<Float> / 2.0f + dr::Pi<Float>*i/N);
            // ScalarFloat arg = a * (dr::Pi<Float>*(i/(2.0f*N) - offsets[angle/5]));
            ScalarFloat arg = a * (dr::Pi<Float>*i/(2.0f*N));
            // ScalarFloat arg = a * (20.0f*i / N);
            ScalarFloat sincvalue = dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo);
            // Mask mask = (dr::isfinite(sincvalue) | theta > offset);
            // Mask mask = dr::isfinite(sincvalue);
            // argが0の時はゼロ除算を防ぐため
            // if ((i / (2.0f*N)) <= offsets[angle/5]) {sincvalue = ScalarFloat(1.0);}
            sincvalue = dr::select(
                dr::isfinite(sincvalue),
                dr::abs(sincvalue),
                ScalarFloat(1.0)
            );
            data[i] = sincvalue;
        }
        // printf("%f\n", data[2]);

        // 正規化（ごり押し）
        float sum = 0.0f;
        for (int i = 0; i < N; ++i)
        {
            sum += data[i];
        }
        // std::cout << "[" << std::endl;
        // float sum = dr::Pi<Float>;
        for (int i = 0; i < N; ++i)
        {
            data[i] /= sum;
            // data[i] *= 0.4;
        }

        sum_rawPDF = sum;
        m_data = &data[0]; // LUTの先頭ポインタ
        struct DiscreteDistribution<Float> dd(m_data, N);
        m_pdfdata = dd;
    }

    void make_sincarray1_multi(float a, int expo)
    {
        ScalarFloat data[N];
        int i = 0, idx = 0;

        for (int i = 0; i < N; ++i)
        {
            ScalarFloat sincvalue = ScalarFloat(0.0f);
            for (int n = -100; n <= 100; ++n) {
                ScalarFloat arg = a * ((dr::Pi<Float>*i / (2.0f*N)) - n*0.08);
                sincvalue += dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo);
                // argが0の時はゼロ除算を防ぐため
                sincvalue = dr::select(
                    dr::isfinite(sincvalue),
                    dr::abs(sincvalue),
                    ScalarFloat(1.0)
                );
            }
            
            data[i] = sincvalue;
        }

        // 正規化（ごり押し）
        float sum = 0.0f;
        for (int i = 0; i < N; ++i)
        {
            sum += data[i];
        }
        // float sum = dr::Pi<Float>;
        for (int i = 0; i < N; ++i)
        {
            data[i] /= sum;
            // data[i] *= 0.4;
        }

        sum_rawPDF = sum;
        m_data = &data[0]; // LUTの先頭ポインタ
        struct DiscreteDistribution<Float> dd(m_data, N);
        m_pdfdata = dd;
    }

    void traverse(TraversalCallback *callback) override
    {
        callback->put_object("divergence", m_divergence.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("angle", m_angle.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("cornersize", m_cornersize.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);

        isTraversed = true;

        // mi.traverse(scene)が呼び出されたときに新しいaのLUTを作る
        // M = 1440; N = 360;
        // make_sincarray2(m_a.get()->max(), 2);

        float a = m_a.get()->max();
        int angle = static_cast<int>(m_angle.get()->max());
        // offset = offsets[angle/5]*dr::Pi<Float>;
        // offset = 0.f;
        // offset = 0.000400*dr::Pi<Float>;
        // float upper = 0.0001555555555554644, lower = 0.000700;
        // float upper = 0.0004000, lower = 0.000700;
        // float upper = 0.001000, lower = 0.001200;
        // float upper = 0.0001555, lower = 0.000255;
        // float upper = 0.0, lower = 0.0;
        // offset = (((upper - lower)/340)*(a - 540) + upper)*dr::Pi<Float>;
        offset = 0.0001555555555554644*dr::Pi<Float>;
        // offset = 0.000020000000006348273*dr::Pi<Float>;
        // offset = dr::select(
        //     offset >= 0.f,
        //     offset,
        //     Float(0.0)
        // );
        // printf("%f\n", offsets[angle/5]);
        make_sincarray1(m_a.get()->max(), 2, angle);
        // make_sincarray1_multi(m_a.get()->max(), 2);

        // offset = 0.f;
        // offset2 = (a+2)*dr::Pi<Float>/(a*(a+1));
        // offset2 = 2*dr::Pi<Float> / a;
        // offset2 = dr::Pi<Float> / a;
        offset2 = 0.f;
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
        // SurfaceInteraction3f &si,
        Float sample1,
        const Point2f &sample2,
        Mask active
        ) const override
        // ) const
    {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float sin_theta_i = Frame3f::sin_theta(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float sin_phi_i = Frame3f::sin_phi(si.wi);
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};

        // Point2u point2;
        // Point2f point2;
        Float point;
        Float pdfvalue;
        Point2f sampled;
        // std::tie(point2, pdfvalue, sampled) = m_pdfdata2.sample(sample2, active);
        // point = m_pdfdata.sample(sample1, active);
        point = m_pdfdata.sample(sample2.y());

        // ----- 出射光のサンプリング -----
        // theta, phiのサンプリング
        // index_theta = point;

        Float del_phi = dr::TwoPi<Float> * sample2.x(); // 角度の差分にする
        // printf("phi = %f\n", (del_phi*180.f)/dr::Pi<Float>);
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N) + dr::Pi<Float>/15555.0f;
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N) + dr::Pi<Float>/(2.0*m_a->eval(si, active).x());
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N) + dr::Pi<Float>*0.0276583f;
        Float a = m_a->eval(si, active).x();
        // Float del_theta = -dr::Pi<Float> / 2.0f + point*dr::Pi<Float>/N + offset + offset2;
        auto rng = setRandomGenerator(sample1*10000000);
        Float random1 = rand(rng);
        // Float del_theta = point*dr::Pi<Float>/(2.0*N) + (offset + offset2)*random1;
        Float del_theta = point*dr::Pi<Float>/(2.0*N) + (offset + offset2);
        // Float del_theta = point*dr::Pi<Float>/(2.0*N) + offset2*sample1;
        // Float tan_theta = point*dr::Pi<Float>/(2.0f*N);
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N);
        // Float del_theta = dr::atan(tan_theta) + (offset + offset2)*sample1;
        // Float del_theta = point*dr::Pi<Float>/(2.0f*N);
        // printf("%f\n", point);
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_theta);
        // Vector3f delvec = Vector3f(0.f, 0.f, 1.f);
        // Vector3f norm = si.n;
        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        // printf("%f, %f, %f\n", si.wi.x(), si.wi.y(), si.wi.z());
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;

        bs.sampled_component = 0;
        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        bs.eta = 1.f;
        bs.pdf = 1.f;

        // BSSRDF
        // auto rng = setRandomGenerator(sample1*1000000);
        Point3f p = si.p;
        Float r1 = rand(rng);
        Float r2 = rand(rng);
        float cornersize = m_cornersize.get()->max();
        Float shift = cornersize*r1;
        Float shift2 = cornersize*r2;
        // Float shift = cornersize;
        Float r3 = rand(rng);
        Float theta = 2*r3*dr::Pi<Float>;
        // Float theta = r3*360.0;
        // Float theta = del_phi + dr::Pi<Float>;
        // Float theta = del_phi;
        // Float theta = del_phi + dr::Pi<Float>;

        // 何故この書き方をしているのか。
        // bs.p = Point3f(shift*dr::cos(theta), shift*dr::sin(theta), 0.0f);
        // bs.p = Point3f(cornersize, 0.0f, 0.0f);
        bs.p = Point3f(cornersize, cornersize, 0.0f);
        // bs.p = Point3f(0.0f, 0.0f, -cornersize);
        // bs.p = Point3f(shift*dr::cos(theta), 0.f, shift*dr::sin(theta));
        // printf("si.p = %f, %f, %f\n", si.p.x(), si.p.y(), si.p.z());
        // bs.p = Point3f(shift*dr::cos(theta), shift2*dr::sin(theta), 0.0);
        // bs.isBSSRDF = true;
        // bs.p = 0.0*si.wi + Point3f(shift*dr::cos(theta), shift*dr::sin(theta), -1.0);
        // bs.p = si.wi;
        // bs.p = p;


        UnpolarizedSpectrum arg = m_a->eval(si, active)*(del_theta - offset - offset2);
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask mask = arg == 0;
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue /= sum_rawPDF;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active);

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

        Float a = m_a->eval(si, active).x();
        UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset - offset2);
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset2);
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask mask = arg == 0;
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue /= sum_rawPDF;
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
        Float delta_theta = acos(cos_delta_theta);

        Float a = m_a->eval(si, active).x();
        // Float offset2 = 2*dr::Pi<Float>/(a*(a+1));
        UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset - offset2);
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset2);
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
        // UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        // Mask mask = arg == 0;
        sincvalue = dr::select(
            dr::isfinite(sincvalue),
            sincvalue,
            ScalarFloat(1.0)
        );
        // sincvalue /= sum_rawPDF;

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue.x(), 0.f);
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

            Float a = m_a->eval(si, active).x();
            // Float offset = dr::Pi<Float>/15555.0f;
            // Float offset2 = 2*dr::Pi<Float>/(a*(a+1));
            UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset - offset2);
            // UnpolarizedSpectrum arg = m_a->eval(si, active)*(delta_theta - offset2);
            UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
            // UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
            // UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
            // Mask mask = arg == 0;
            sincvalue = dr::select(
                dr::isfinite(sincvalue),
                sincvalue,
                ScalarFloat(1.0)
            );
            // sincvalue /= sum_rawPDF;
            UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue;
            Float pdf = dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue.x(), 0.f);

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
        }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "sincarray2[" << std::endl
            << " divergence = " << string::indent(m_divergence) << std::endl
            << " reflectance = " << string::indent(m_reflectance) << std::endl
            << "a = " << string::indent(m_a) << std::endl
            << "]";
        return oss.str();
    }

    // std::string to_string() const override {
    //     std::ostringstream oss;
    //     oss << "a = " << m_a.get() << " .json" << std::endl;
    //     cout << "a = " << m_a.get() << std::endl;
    //     return oss.str();
    // }

    std::string make_LUTPath() const {
        std::ostringstream oss;
        std::ostringstream out;
        int precision = 2;
        out << std::fixed << std::setprecision(precision) << m_a.get()->max();
        std::string a_str = out.str();
        oss << LUTDir << "pdf_sinc2, a = " << a_str << ".json";
        return oss.str();
    }


    MI_DECLARE_CLASS()

// フィールド
private:
    uint32_t M; // LUTのphiサイズ
    uint32_t N; // LUTのthetaサイズ
    float sum_rawPDF;
    Point2u point2;
    ref<Texture> m_divergence;
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    vector<vector<float>> m_beforedata;
    ScalarFloat *m_data;
    DiscreteDistribution<Float> m_pdfdata;
    DiscreteDistribution2D<Float, 2> m_pdfdata2;
    // Hierarchical2D<Float, 2> m_pdfdata;
    ScalarVector2u m_size;
    std::string LUTFilename;
    std::string LUTDir;
    Float index_phi, index_theta;
    bool isTraversed;
    Float offset;
    Float offset2;
    float offsets[10] = {0.0007000f, 0.0001555555555554644f, 0.0001555555555554644f, 0.0001555555555554644f, 0.0001555555555554644f, 0.0001555555555554644f, 0.0001555555555554644f, 0.0001555555555554644f, 0.0001555555555554644f, 0.0001200f};
    // float offsets[10] = {0.004023861111, 0.004023861111, 0.004023861111, 0.004023861111, 0.004023861111, 0.004023861111, 0.004023861111, 0.004023861111, 0.004023861111, 0.004023861111};
    ref<Texture> m_angle;
    ref<Texture> m_cornersize;
};

MI_IMPLEMENT_CLASS_VARIANT(sincarray2, BSDF)
MI_EXPORT_PLUGIN(sincarray2, "sincarray2")
NAMESPACE_END(mitsuba)