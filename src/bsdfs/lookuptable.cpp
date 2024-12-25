#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <random>
#include "nlohmann/json.hpp"

NAMESPACE_BEGIN(mitsuba)

using namespace std;
using json = nlohmann::json;

/**!

.. _bsdf-diffuse:

Smooth diffuse material (:monosp:`diffuse`)
-------------------------------------------

.. pluginparameters::

 * - reflectance
   - |spectrum| or |texture|
   - Specifies the diffuse albedo of the material (Default: 0.5)
   - |exposed|, |differentiable|

The smooth diffuse material (also referred to as *Lambertian*)
represents an ideally diffuse material with a user-specified amount of
reflectance. Any received illumination is scattered so that the surface
looks the same independently of the direction of observation.

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_plain.jpg
   :caption: Homogeneous reflectance
.. subfigure:: ../../resources/data/docs/images/render/bsdf_diffuse_textured.jpg
   :caption: Textured reflectance
.. subfigend::
   :label: fig-diffuse

Apart from a homogeneous reflectance value, the plugin can also accept
a nested or referenced texture map to be used as the source of reflectance
information, which is then mapped onto the shape based on its UV
parameterization. When no parameters are specified, the model uses the default
of 50% reflectance.

Note that this material is one-sided---that is, observed from the
back side, it will be completely black. If this is undesirable,
consider using the :ref:`twosided <bsdf-twosided>` BRDF adapter plugin.
The following XML snippet describes a diffuse material,
whose reflectance is specified as an sRGB color:

.. tabs::
    .. code-tab:: xml
        :name: diffuse-srgb

        <bsdf type="diffuse">
            <rgb name="reflectance" value="0.2, 0.25, 0.7"/>
        </bsdf>

    .. code-tab:: python

        'type': 'diffuse',
        'reflectance': {
            'type': 'rgb',
            'value': [0.2, 0.25, 0.7]
        }

Alternatively, the reflectance can be textured:

.. tabs::
    .. code-tab:: xml
        :name: diffuse-texture

        <bsdf type="diffuse">
            <texture type="bitmap" name="reflectance">
                <string name="filename" value="wood.jpg"/>
            </texture>
        </bsdf>

    .. code-tab:: python

        'type': 'diffuse',
        'reflectance': {
            'type': 'bitmap',
            'filename': 'wood.jpg'
        }
*/
template <typename Float, typename Spectrum>
class LookupTable final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    LookupTable(const Properties &props) : Base(props) {
        m_reflectance = props.texture<Texture>("reflectance", .5f);
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);

        // LUTを読み込む
        // vector<vector<float>> data;
        std::string dirpath = "/home/sugawara.ryo/BRDFEstimation/json/";
        // std::string filename = dirpath+"pdf_LinearDist.json";
        // std::string filename_cum = dirpath+"cum_LinearDist.json";
        // std::string filename_sample = dirpath+"samples_LinearDist.json";
        std::string filename = dirpath+"pdf_peakDist.json";
        std::string filename_cum = dirpath+"cum_peakDist.json";
        std::string filename_sample = dirpath+"samples_peakDist.json";
        // std::string filename = dirpath+"pdf_bokeDist.json";
        // std::string filename_cum = dirpath+"cum_bokeDist.json";
        // std::string filename_sample = dirpath+"samples_bokeDist.json";
        ifstream ifs(filename.c_str());
        if (ifs.good())
        {
            json m_json;
            ifs >> m_json;

            for (const auto& json : m_json)
            {
                vector<float> d = json.get<vector<float>>();
                data.push_back(d);
            }
        }
        else
        {
            cout << "ファイルの読み込みに失敗しました" << endl;
        }

        // 周辺化を行う
        vector<float> marginalized_y = marginalize(data, 1);
        // yの累積分布関数を求める
        cum_phi = cumulative_sum(marginalized_y);
        // xごとの累積分布関数を求める
        
        cout << "json is read" << endl;
        cout << "test : " << data[1][5] << endl;
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("reflectance", m_reflectance.get(), +ParamFlags::Differentiable);
    }

    // 確率密度関数を積分したら1になるように正規化する
    vector<vector<float>>& normalize(vector<vector<float>>& fvector)
    {
        float sum = 0.0f;
        for (auto& parts : fvector)
        {
            float tempsum = accumulate(parts.begin(), parts.end(), remove_reference<decltype(*std::begin(parts))>::type{0});
            sum += tempsum;
        }
        std::transform(cbegin(fvector), cend(fvector), begin(fvector), [&sum](vector<float> vec) {
            std::transform(cbegin(vec), cend(vec), begin(vec), [&sum](float v) {
                return v / sum;
            });
            return vec;
        });

        return fvector;
    }


    vector<float>& normalize(vector<float>& fvector) const
    {
        float sum = 0.0f;
        for (auto& val : fvector)
        {
            sum += val;
        }
        std::transform(cbegin(fvector), cend(fvector), begin(fvector), [&sum](float v) {
            return v / sum;
        });

        return fvector;
    }

    // テンプレート型にした方が良いかもしれない
    // 累積分布関数を作る
    vector<float> cumulative_sum(vector<float>& pdf) const
    {
        vector<float> result(pdf.size());
        partial_sum(pdf.begin(), pdf.end(), result.begin());
        return result;
    }

    // 周辺化する
    // axis ... 0のときはxの関数として返す、1のときはyの関数として返す
    vector<float> marginalize(vector<vector<float>>& pdf, int axis)
    {
        vector<float> result(pdf.size());

        switch (axis)
        {
            case 0:
                // yで積分
                break;
            case 1:
                // xで積分
                {
                    int i = 0;
                    for (auto& xvec : pdf)
                    {
                        float temp = accumulate(xvec.begin(), xvec.end(), remove_reference<decltype(*std::begin(xvec))>::type{});
                        result[i++] = temp;
                    }
                    break;
                }
            default:
                break;
        }

        return result;
    }

    // [0, 1]の範囲で乱数を生成
    float generate_random() const
    {
        random_device seed_gen;
        static std::mt19937 generator(seed_gen());
        static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
        return distribution(generator);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* sample1 */,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f; // 入射光が表面上にあるか
        // ativeでない場合は0として返す
        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
            return { bs, 0.f };


        // ----- 出射光のサンプリング -----
        // 最初にphiのサンプリング
        float u_phi = generate_random();
        auto itr_u_phi = lower_bound(cum_phi.begin(), cum_phi.end(), u_phi);
        auto index_u_phi = distance(cum_phi.begin(), itr_u_phi);
        // つづいてthetaのサンプリング
        vector<float> pdf_theta = data[index_u_phi]; // thetaの分布を取得
        vector<float> pdf_thetanorm = normalize(pdf_theta); // 正規化
        vector<float> cum_x = cumulative_sum(pdf_thetanorm);
        float v_theta = generate_random();
        auto itr_v_theta = lower_bound(cum_x.begin(), cum_x.end(), v_theta);
        auto index_v_theta = distance(cum_x.begin(), itr_v_theta);
        vector<Float> sampled_angles{index_u_phi*dr::TwoPi<Float>/cum_phi.size(), index_v_theta*dr::Pi<Float>/(2.0f*pdf_theta.size())}; // サンプリングしたphi, thetaのペア
        sampled_angles[0] -= dr::Pi<Float>; sampled_angles[1] -= dr::Pi<Float>/2.0f; // 角度の差分にする

        cout << sampled_angles << endl;
        Float wo_x = dr::cos(sampled_angles[0])*dr::sin(sampled_angles[1]);
        Float wo_y = dr::sin(sampled_angles[0])*dr::sin(sampled_angles[1]);
        Float wo_z = dr::cos(sampled_angles[1]);
        // Normal3f h = Vector3f(wo_x, wo_y, wo_z);

        // bs.wo = warp::square_to_cosine_hemisphere(sample2); // 出射光のサンプリング
        Vector3f ideal_wo = reflect(si.wi);
        bs.wo = Vector3f(wo_x, wo_y, wo_z) + ideal_wo;
        // bs.wo  = reflect(si.wi, h); // Rusinkiewicz
        // bs.wo = ideal_wo;
        // bs.pdf = warp::square_to_cosine_hemisphere_pdf(bs.wo); // 確率密度関数の値
        bs.pdf = data[index_u_phi][index_v_theta];
        bs.eta = 1.f; 
        bs.sampled_type = +BSDFFlags::DiffuseReflection;
        bs.sampled_component = 0;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);

        return { bs, depolarizer<Spectrum>(value) & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f; // 表面化に行ったら無効

        // レンダリング方程式の値
        UnpolarizedSpectrum value =
            m_reflectance->eval(si, active) * dr::InvPi<Float> * cos_theta_o;

        return depolarizer<Spectrum>(value) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        // Float pdf = warp::square_to_cosine_hemisphere_pdf(wo);
        Float pdf = data[index_u_phi][index_v_theta];

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f); // ３項演算子と同じ役割
    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return { 0.f, 0.f };

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        UnpolarizedSpectrum value =
            m_reflectance->eval(si, active) * dr::InvPi<Float> * cos_theta_o;

        Float pdf = data[index_u_phi][index_v_theta];

        return { depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f) };
    }

    Spectrum eval_diffuse_reflectance(const SurfaceInteraction3f &si,
                                      Mask active) const override {
        return m_reflectance->eval(si, active);
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "LookupTable[" << std::endl
            << "  reflectance = " << string::indent(m_reflectance) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
    vector<vector<float>> data;
    vector<float> cum_phi;
    int index_u_phi, index_v_theta;
};

MI_IMPLEMENT_CLASS_VARIANT(LookupTable, BSDF)
MI_EXPORT_PLUGIN(LookupTable, "Smooth diffuse material")
NAMESPACE_END(mitsuba)
