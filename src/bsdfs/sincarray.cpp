#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/texture.h>

#include <drjit/dynamic.h>
#include <mitsuba/core/distr_2d.h>

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
class sincarray final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    sincarray(const Properties &props) : Base(props)
    {
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        m_divergence = props.texture<Texture>("divergence", 1.f);
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        if (props.has_property("LUTFilename"))
        {
            LUTFilename = props.string("LUTFilename");
        }
        else
        {
            Throw("The property 'filename' is required!");
        }

        if (props.has_property("LUTDir"))
        {
            LUTDir = props.string("LUTDir");
        }
        else
        {
            Throw("The property 'LUTDir' is required!");
        }

        m_a = props.texture<Texture>("a", 1.f);

        isTraversed = false;

        // std::string Filename = LUTDir + LUTFilename;
        // std::cout << Filename << std::endl;
        // ifstream ifs(Filename);
        // if (ifs.good())
        // {
        //     json m_json;
        //     ifs >> m_json;

        //     int j = 0;
        //     for (const auto& json : m_json)
        //     {
        //         vector<float> d_float = json.get<vector<float>>();
        //         m_beforedata.push_back(d_float);
        //     }
        //     M = m_beforedata[0].size();
        //     N = m_beforedata.size();

        //     ScalarVector2u size(M, N);
        //     ScalarFloat data[M*N];
        //     int i = 0, idx = 0;
        //     for (const auto& l : m_beforedata)
        //     {
        //         idx = M*i;
        //         for (const auto& v : l)
        //         {
        //             data[idx++] = v*ScalarFloat(1.0);
        //         }
        //         ++i;
        //     }

        //     m_data = &data[0];
        //     m_size = size;
        //     m_pdfdata = DiscreteDistribution2D<Float, 2>(&data[0], size);

        //     std:cout << "Initialized：ファイルを読んだ" << Filename << std::endl;
        // }
        // else
        // {
        //     std::cout << "ファイルの読み込みに失敗しました" << std::endl;
        // }
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

    vector<vector<float>> make_sincarray(vector<vector<float>>& vecs, int M, int N, float a, int expo) const
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

    void traverse(TraversalCallback *callback) override
    {
        LUTFilename = make_LUTPath();
        std::cout << "LUTFilename = " << LUTFilename << std::endl;

        callback->put_object("divergence", m_divergence.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);

        isTraversed = true;

        // mi.traverse(scene)が呼び出されたときに新しいaのLUTを作る
        if (isTraversed)
        {
            // もともとあるlUTを一旦削除
            for (auto& l : m_beforedata)
            {
                l.clear();
            }
            m_beforedata.clear();
            // 新たなLUTを読み込む
            std::cout << LUTFilename << std::endl;
            ifstream ifs(LUTFilename);
            if (ifs.good())
            {
                json m_json;
                ifs >> m_json;

                int j = 0;
                for (const auto& json : m_json)
                {
                    vector<float> d_float = json.get<vector<float>>();
                    m_beforedata.push_back(d_float);
                }
                M = m_beforedata[0].size();
                N = m_beforedata.size();

                ScalarVector2u size(M, N);
                ScalarFloat data[M*N];
                int i = 0, idx = 0;
                for (const auto& l : m_beforedata)
                {
                    idx = M*i;
                    for (const auto& v : l)
                    {
                        data[idx++] = v*ScalarFloat(1.0);
                    }
                    ++i;
                }

                m_data = &data[0];
                m_size = size;
                m_pdfdata = DiscreteDistribution2D<Float, 2>(&data[0], size);

                std:cout << "Traverse：ファイルを読んだ" << LUTFilename << std::endl;
            }
            else
            {
                std::cout << "ファイルの読み込みに失敗しました" << std::endl;
            }
        }
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
        active &= cos_theta_i > 0.f;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};

        Point2u point2;
        Float pdfvalue;
        Point2f sampled;
        std::tie(point2, pdfvalue, sampled) = m_pdfdata.sample(sample2, active);

        // ----- 出射光のサンプリング -----
        // 最初にthetaのサンプリング
        auto index_theta = point2.y();
        // つづいてthetaのサンプリング
        auto index_phi = point2.x();
        // 修正後
        Float sampled_angles[2] = {point2.x()*dr::TwoPi<Float>/M, dr::abs(point2.y()*dr::Pi<Float>/N - dr::Pi<Float>/2.0f)}; // サンプリングしたphi, thetaのペア
        sampled_angles[0] = dr::TwoPi<Float> * sample1; // 角度の差分にする

        Float del_phi = sampled_angles[0], del_theta = sampled_angles[1];
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_phi);
        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;
        // bs.wo = si.wi;

        bs.sampled_component = 0;
        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        bs.eta = 1.f;
        bs.pdf = 1.f;

        UnpolarizedSpectrum arg = m_a->eval(si, active)*del_theta;
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            1.0f,
            sincvalue
        );

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
        Float cos_phi_i = Frame3f::cos_phi(si.wi),
              cos_phi_o = Frame3f::cos_phi(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        Float theta_i = dr::acos(cos_theta_i);
        Float theta_o = dr::acos(cos_theta_o);
        Float delta_theta = theta_o - theta_i;
        
        Float phi_i = dr::acos(cos_phi_i);
        Float phi_o = dr::acos(cos_phi_o);
        Float delta_phi = phi_o - phi_i;

        UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            1.0f,
            sincvalue
        );
        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active);

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
        Float cos_phi_i = Frame3f::cos_phi(si.wi),
              cos_phi_o = Frame3f::cos_phi(wo);
        
        Float theta_i = dr::acos(cos_theta_i);
        Float theta_o = dr::acos(cos_theta_o);
        Float delta_theta = theta_o - theta_i;
        
        Float phi_i = dr::acos(cos_phi_i);
        Float phi_o = dr::acos(cos_phi_o);
        Float delta_phi = phi_o - phi_i;

        UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            1.0f,
            sincvalue
        );

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue[0], 0.f);
        // return 1.0f;
        // return 0.f;
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
            Float cos_phi_i = Frame3f::cos_phi(si.wi),
                  cos_phi_o = Frame3f::cos_phi(wo);

            Float theta_i = dr::acos(cos_theta_i);
            Float theta_o = dr::acos(cos_theta_o);
            Float delta_theta = theta_o - theta_i;
            
            Float phi_i = dr::acos(cos_phi_i);
            Float phi_o = dr::acos(cos_phi_o);
            Float delta_phi = phi_o - phi_i;

            active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

            UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
            UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
            sincvalue = dr::select(
                dr::isnan(sincvalue),
                1.0f,
                sincvalue
            );
            // UnpolarizedSpectrum value = m_reflectance->eval(si, active);
            UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue;
            Float pdf = dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, sincvalue[0], 0.f);
            // Float pdf = sincvalue[0];

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
        }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "sincarray[" << std::endl
            << " divergence = " << string::indent(m_divergence) << std::endl
            << " reflectance = " << string::indent(m_reflectance) << std::endl
            << " LUTFilename = " << LUTFilename << std::endl
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
    Point2u point2;
    ref<Texture> m_divergence;
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    vector<vector<float>> m_beforedata;
    ScalarFloat *m_data;
    DiscreteDistribution2D<Float, 2> m_pdfdata;
    ScalarVector2u m_size;
    std::string LUTFilename;
    std::string LUTDir;
    int index_phi, index_theta;
    bool isTraversed;
};

MI_IMPLEMENT_CLASS_VARIANT(sincarray, BSDF)
MI_EXPORT_PLUGIN(sincarray, "sincarray")
NAMESPACE_END(mitsuba)