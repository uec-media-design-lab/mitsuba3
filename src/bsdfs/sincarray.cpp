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

// TODOリスト
// ・FloatのLookup Tableを作るなら、drjitに定義されている3次元配列（Floatの2次元配列）を使うようにする？外部公開する時向けか

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class sincarray final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    // using ScalarFloat = dr::scalar_t<Float>;
    // using ScalarFloat = Float;
    // using FloatStorage = DynamicBuffer<Float>;
    // Phase Function
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

        LUTFilename = "/home/sugawara.ryo/BRDFEstimation/pdf_sinc2, a = 10.00.json";
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
            // ScalarVector2u size(N);
            // Float data[M*N];
            // Float *data;
            // ScalarFloat *data;
            ScalarFloat data[M*N];
            // std::vector<Float> data;
            int i = 0, idx = 0;
            for (const auto& l : m_beforedata)
            {
                idx = M*i;
                for (const auto& v : l)
                {
                    data[idx++] = v*ScalarFloat(1.0);
                    // data[idx++] = v*Float(1.0);
                    // data[idx++] = v;
                    // data.push_back(v);
                }
                ++i;
                // std::cout << "i = " << i << ", idx = " << idx << std::endl;
            }

            m_data = &data[0];
            m_size = size;
            // m_pdfdata = DiscreteDistribution2D(data, size);
            m_pdfdata = DiscreteDistribution2D<Float, 2>(&data[0], size);
            // m_pdfdata = DiscreteDistribution2D<Float, 2>(data, size);


            std:cout << "ファイルを読んだ" << LUTFilename << std::endl;
        }
        else
        {
            std::cout << "ファイルの読み込みに失敗しました" << std::endl;
        }

        // 周辺化を行う
        // vector<float> marginalized_y = marginalize(m_beforedata, 1);
        // // vector<Float> marginalized_y = marginalize(m_beforedata, 1);
        // // yの累積分布関数を求める
        // cum_theta = cumulative_sum(marginalized_y);
        // xごとの累積分布関数を求める
        
        std::cout << "json is read" << std::endl;
        std::cout << "test : " << m_data[5] << std::endl;
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

        // // x軸周りの回転
        // Float temp_y = sin_theta_i*vec.y();
        // Float temp_z = cos_theta_i*vec.z();

        // // z軸周りの回転
        // Float result_x = cos_phi_i*vec.x();
        // Float result_y = sin_phi_i*temp_y;
        // Float result_z = temp_z;
        // ロドリゲスの回転公式
        Vector3f fromNorm = dr::normalize(normal);
        Vector3f toNorm = dr::normalize(to);
        Float s = mag(dr::cross(fromNorm, toNorm));
        Float c = dr::dot(fromNorm, toNorm);
        Vector3f result = c*vec + dr::dot(axis, vec)*(1-c)*axis + dr::cross(axis, vec)*s;
        result = dr::normalize(result);
        // cout << result << endl;
        // return Vector3f(result_x, result_y, result_z);

        return result;
    }

    vector<vector<float>> make_sincarray(vector<vector<float>>& vecs, int M, int N, float a, int expo) const
    {
        // vector<vector<float>> v;
        for (int i = 0; i < N; i++) { 
            vector<float> parts = vector<float>(M); // 0で初期化
            for (int j = 0; j < M; j++)
            {
                float arg = a * static_cast<float>(std::sqrt( std::pow(j - M/2, 2) + std::pow(i - N/2, 2) ));
                parts[j] = arg == 0 ? 1.0f : std::pow(std::sin(arg), expo) / std::pow(arg, expo); // argが0の時はゼロ除算を防ぐため
                cout << "parts["<< j << "] = " << parts[j] << endl;
            }
            vecs.push_back(parts);
            // v[i] = parts;
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
        // callback->put_object("a", m_a.get(), +ParamFlags::Differentiable); // discontinuousを消した

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
                // ScalarVector2u size(N);
                // Float data[M*N];
                // Float *data;
                // ScalarFloat *data;
                ScalarFloat data[M*N];
                // std::vector<Float> data;
                int i = 0, idx = 0;
                for (const auto& l : m_beforedata)
                {
                    idx = M*i;
                    for (const auto& v : l)
                    {
                        data[idx++] = v*ScalarFloat(1.0);
                        // data[idx++] = v*Float(1.0);
                        // data[idx++] = v;
                        // data.push_back(v);
                    }
                    ++i;
                    // std::cout << "i = " << i << ", idx = " << idx << std::endl;
                }

                m_data = &data[0];
                m_size = size;
                // m_pdfdata = DiscreteDistribution2D(data, size);
                m_pdfdata = DiscreteDistribution2D<Float, 2>(&data[0], size);
                // m_pdfdata = DiscreteDistribution2D<Float, 2>(data, size);


                std:cout << "ファイルを読んだ" << LUTFilename << std::endl;
            }
            else
            {
                std::cout << "ファイルの読み込みに失敗しました" << std::endl;
            }

            // 周辺化を行う
            // vector<float> marginalized_y = marginalize(m_beforedata, 1);
            // // vector<Float> marginalized_y = marginalize(m_beforedata, 1);
            // // yの累積分布関数を求める
            // cum_theta = cumulative_sum(marginalized_y);
            // xごとの累積分布関数を求める
            
            std::cout << "json is read" << std::endl;
            std::cout << "test : " << m_data[5] << std::endl;
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

        // std::cout << "a = " << m_a << std::endl; // TODO消す

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        Float sin_theta_i = Frame3f::sin_theta(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float sin_phi_i = Frame3f::sin_phi(si.wi);
        active &= cos_theta_i > 0.f;

        // std::cout << "cos_theta_i = " << cos_theta_i << std::endl;

        BSDFSample3f bs = dr::zeros<BSDFSample3f>();
        if (unlikely(dr::none_or<false>(active) || !ctx.is_enabled(BSDFFlags::DiffuseReflection))) return {bs, .0f};

        Point2u point2;
        Float pdfvalue;
        Point2f sampled;
        std::tie(point2, pdfvalue, sampled) = m_pdfdata.sample(sample2, active);

        // ----- 出射光のサンプリング -----
        // 最初にthetaのサンプリング
        // auto itr_v_theta = lower_bound(cum_theta.begin(), cum_theta.end(), generate_random());
        auto index_theta = point2.y();
        // // つづいてthetaのサンプリング　書き換え前
        // // 書き換え後
        // auto pdf_phi = m_beforedata[index_theta]; // phiの分布を取得
        // auto pdf_phinorm = normalize(pdf_phi); // 正規化
        // auto cum_phi = cumulative_sum(pdf_phinorm);
        // auto itr_u_phi = lower_bound(cum_phi.begin(), cum_phi.end(), generate_random());
        auto index_phi = point2.x();
        // Float sampled_angles[2] = {point2.x()*dr::TwoPi<Float>/M, point2.y()*dr::Pi<Float>/(2.0f*N)}; // サンプリングしたphi, thetaのペア
        // sampled_angles[0] = dr::TwoPi<Float> * generate_random(); // 角度の差分にする
        // // sampled_angles[0] = dr::TwoPi<Float> * sample1; // 角度の差分にする
        // sampled_angles[1] -= dr::Pi<Float>/4.0f; // 角度の差分にする
        // 修正後
        Float sampled_angles[2] = {point2.x()*dr::TwoPi<Float>/M, dr::abs(point2.y()*dr::Pi<Float>/N - dr::Pi<Float>/2.0f)}; // サンプリングしたphi, thetaのペア
        sampled_angles[0] = dr::TwoPi<Float> * sample1; // 角度の差分にする

        // cout << "phi = " << sampled_angles[0] << ", theta = " << sampled_angles[1] << endl;
        Float del_phi = sampled_angles[0], del_theta = sampled_angles[1];
        // std::cout << "del_phi = " << del_phi << ", del_theta = " << del_theta << std::endl;
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_phi);
        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;
        // bs.wo = si.wi;

        bs.sampled_component = 0;
        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        // bs.wo = 
        bs.eta = 1.f;
        bs.pdf = 1.f;

        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * m_pdfdata.eval(point2, active);
        // UnpolarizedSpectrum arg = m_a->eval(si, active)*dr::sqrt(dr::pow(del_phi, 2) + dr::pow(del_theta, 2));
        UnpolarizedSpectrum arg = m_a->eval(si, active)*del_theta;
        // UnpolarizedSpectrum sincvalue = dr::select(arg == Float(0), dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2), Float(1.0)); // 0除算を防ぐため
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            1.0f,
            sincvalue
        );
        // sincvalue = Float(1.0);

        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * pdfvalue;
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active);
        // UnpolarizedSpectrum value = sincvalue;

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

        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * m_beforedata[index_phi][index_theta];
        // UnpolarizedSpectrum value = m_reflectance->eval(si, active) * m_pdfdata.eval(point2, active);

        Float theta_i = dr::acos(cos_theta_i);
        Float theta_o = dr::acos(cos_theta_o);
        Float delta_theta = theta_o - theta_i;
        
        Float phi_i = dr::acos(cos_phi_i);
        Float phi_o = dr::acos(cos_phi_o);
        Float delta_phi = phi_o - phi_i;

        // UnpolarizedSpectrum arg = m_a->eval(si, active)*dr::sqrt(dr::pow(delta_phi, 2) + dr::pow(delta_theta, 2));
        UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
        // UnpolarizedSpectrum sincvalue = dr::select(dr::any(arg != 0.f), dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2), Float(1.0));
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            1.0f,
            sincvalue
        );
        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * sincvalue;

        // UnpolarizedSpectrum value = m_reflectance->eval(si, active);
        // UnpolarizedSpectrum value = sincvalue;
        return depolarizer<Spectrum>(value) & active;

        // return 0.f;
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

        // UnpolarizedSpectrum arg = m_a->eval(si, active)*dr::sqrt(dr::pow(delta_phi, 2) + dr::pow(delta_theta, 2));
        UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
        // UnpolarizedSpectrum sincvalue = dr::select(dr::any(arg != 0.f), dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2), Float(1.0));
        UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
        sincvalue = dr::select(
            dr::isnan(sincvalue),
            1.0f,
            sincvalue
        );
        // sincvalue = dr::select(!dr::isinf(sincvalue) & !dr::isnan(sincvalue), sincvalue, arg + 1.0);
        // sincvalue = Float(1.0);
        
        // Float pdf = m_beforedata[index_phi][index_theta];
        // Float pdf = m_pdfdata.pdf(point2, active);
        // Float pdf = m_reflectance->eval(si, active);

        // return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
        // return 1.0f;
        // return 0.f;
        return sincvalue[0];
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

            // UnpolarizedSpectrum arg = m_a->eval(si, active)*dr::sqrt(dr::pow(delta_phi, 2) + dr::pow(delta_theta, 2));
            UnpolarizedSpectrum arg = m_a->eval(si, active)*delta_theta;
            // UnpolarizedSpectrum sincvalue = dr::select(dr::any(arg != 0.f), dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2), Float(1.0));
            UnpolarizedSpectrum sincvalue = dr::pow(dr::sin(arg), 2) / dr::pow(arg, 2);
            sincvalue = dr::select(
            dr::isnan(sincvalue),
            1.0f,
            sincvalue
        );
            // sincvalue = dr::select(!dr::isinf(sincvalue) & !dr::isnan(sincvalue), sincvalue, arg + 1.0);
            // sincvalue = Float(1.0);

            // UnpolarizedSpectrum value = m_reflectance->eval(si, active);
            // UnpolarizedSpectrum value = sincvalue;
            UnpolarizedSpectrum value = m_reflectance->eval(si, active)*sincvalue;

            // Float pdf = m_beforedata[index_phi][index_theta];
            // Float pdf = 0.f;
            // Float pdf = m_pdfdata.pdf(point2, active);
            // Float pdf = 1.0f;
            Float pdf = sincvalue[0];

            return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
            // return {0.f, 0.f};
        }

    // std::string to_string() const override {
    //     std::ostringstream oss;
    //     oss << "sincarray[" << std::endl
    //         << " divergence = " << string::indent(m_divergence) << std::endl
    //         << " reflectance = " << string::indent(m_reflectance) << std::endl
    //         << " LUTFilename = " << LUTFilename << std::endl
    //         // << "a = " << m_a.get() << std::endl
    //         << "]";
    //     return oss.str();
    // }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "a = " << m_a.get() << " .json" << std::endl;
        cout << "a = " << m_a.get() << std::endl;
        return oss.str();
    }

    std::string make_LUTPath() const {
        std::ostringstream oss;
        // oss << "/home/sugawara.ryo/BRDFEstimation/json/pdf_sincarray2, a = " << m_a.get() << ".json" << std::endl;
        std::ostringstream out;
        int precision = 2;
        out << std::fixed << std::setprecision(precision) << m_a.get()->max();
        std::string a_str = out.str();
        // oss << "/home/sugawara.ryo/BRDFEstimation/json/pdf_sincarray2, a = " << a_str << ".json";
        oss << LUTDir << "pdf_sinc2, a = " << a_str << ".json";
        return oss.str();
    }


    MI_DECLARE_CLASS()

// フィールド
private:
    uint32_t M; // LUTのphiサイズ
    uint32_t N; // LUTのthetaサイズ
    Point2u point2; // 
    ref<Texture> m_divergence;
    ref<Texture> m_reflectance;
    ref<Texture> m_a;
    vector<vector<float>> m_beforedata;
    // Float *m_data;
    ScalarFloat *m_data;
    DiscreteDistribution2D<Float, 2> m_pdfdata;
    // std::unique_ptr<ScalarFloat[]> m_data;
    // FloatStorage m_data;
    ScalarVector2u m_size;
    // vector<vector<Float>> m_beforedata;
    vector<float> cum_theta;
    std::string LUTFilename;
    std::string LUTDir;
    int index_phi, index_theta;

    bool isTraversed;
};

MI_IMPLEMENT_CLASS_VARIANT(sincarray, BSDF)
MI_EXPORT_PLUGIN(sincarray, "sincarray")
NAMESPACE_END(mitsuba)