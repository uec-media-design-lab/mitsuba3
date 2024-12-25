#include <mitsuba/core/properties.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/ior.h>
#include <mitsuba/render/texture.h>

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
class sinc final : public BSDF<Float, Spectrum>
{
// 関数（コンストラクタやデストラクタ、その他）
public:
    // Phase Function
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    vector<vector<float>> normalize(vector<vector<float>>& fvector) const
    {
        // vector<vector<float>> result;
        float sum = 0.0f;
        // cout << "executing" << endl;
        int i = 0;
        for (const auto& parts : fvector)
        {
            // cout << "size of parts = " << parts.size() << endl;
            // float tempsum = accumulate(parts.begin(), parts.end(), remove_reference<decltype(*std::begin(parts))>::type{0});
            for (const auto& v : parts)
            {
                // cout << "sum = " << sum << ", v = " << v << endl;
                sum += v;
            }
            // float tempsum = std::reduce(parts.begin(), parts.end());
            // cout << "tempsum = " << tempsum << endl;
            // sum += tempsum;
        }
        std::transform(cbegin(fvector), cend(fvector), begin(fvector), [&sum](vector<float> vec) {
            std::transform(cbegin(vec), cend(vec), begin(vec), [&sum](float v) {
                // cout << "v = " << v << ", v/sum = " << v / sum << endl;
                return v / sum;
            });
            return vec;
        });

        // for (const auto& vec : fvector)
        // {
        //     for(const auto& v : vec)
        //     {
        //         cout << "normalizing" << v << endl;
        //     }
        // }

        return fvector;
    }

    vector<float> normalize(vector<float>& fvector) const
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
        // cout << "pdf.size() = " << pdf.size() << endl;
        vector<float> result(pdf.size());
        partial_sum(pdf.begin(), pdf.end(), result.begin());

        // TODO　消す
        // for (const auto& v : result)
        // {
        //     cout << "v of cum_pdf = " << v << endl;
        // }
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
                        float temp = 0.0f;
                        for (const auto& v : xvec)
                        {
                            temp += v;
                            // std::cout << "v = " << v << std::endl;
                            // std::cout << "temp = " << temp << std::endl;
                        }
                        // float temp = accumulate(xvec.begin(), xvec.end(), remove_reference<decltype(*std::begin(xvec))>::type{});
                        result[i++] = temp;
                    }
                    break;
                }
            default:
                break;
        }

        return result;
    }



    // vector<vector<Float>> normalize(vector<vector<Float>>& fvector) const
    // {
    //     // vector<vector<Float>> result;
    //     Float sum = 0.0f;
    //     // cout << "executing" << endl;
    //     int i = 0;
    //     for (const auto& parts : fvector)
    //     {
    //         // cout << "size of parts = " << parts.size() << endl;
    //         // Float tempsum = accumulate(parts.begin(), parts.end(), remove_reference<decltype(*std::begin(parts))>::type{0});
    //         for (const auto& v : parts)
    //         {
    //             // cout << "sum = " << sum << ", v = " << v << endl;
    //             sum += v;
    //         }
    //         // Float tempsum = std::reduce(parts.begin(), parts.end());
    //         // cout << "tempsum = " << tempsum << endl;
    //         // sum += tempsum;
    //     }
    //     std::transform(cbegin(fvector), cend(fvector), begin(fvector), [&sum](vector<Float> vec) {
    //         std::transform(cbegin(vec), cend(vec), begin(vec), [&sum](Float v) {
    //             // cout << "v = " << v << ", v/sum = " << v / sum << endl;
    //             return v / sum;
    //         });
    //         return vec;
    //     });

    //     // for (const auto& vec : fvector)
    //     // {
    //     //     for(const auto& v : vec)
    //     //     {
    //     //         cout << "normalizing" << v << endl;
    //     //     }
    //     // }

    //     return fvector;
    // }

    // vector<Float> normalize(vector<Float>& fvector) const
    // {
    //     Float sum = 0.0f;
    //     for (auto& val : fvector)
    //     {
    //         sum += val;
    //     }
    //     std::transform(cbegin(fvector), cend(fvector), begin(fvector), [&sum](Float v) {
    //         return v / sum;
    //     });

    //     return fvector;
    // }

    // // テンプレート型にした方が良いかもしれない
    // // 累積分布関数を作る
    // vector<Float> cumulative_sum(vector<Float>& pdf) const
    // {
    //     // cout << "pdf.size() = " << pdf.size() << endl;
    //     vector<Float> result(pdf.size());
    //     partial_sum(pdf.begin(), pdf.end(), result.begin());

    //     // TODO　消す
    //     // for (const auto& v : result)
    //     // {
    //     //     cout << "v of cum_pdf = " << v << endl;
    //     // }
    //     return result;
    // }

    // // 周辺化する
    // // axis ... 0のときはxの関数として返す、1のときはyの関数として返す
    // vector<Float> marginalize(vector<vector<Float>>& pdf, int axis)
    // {
    //     vector<Float> result(pdf.size());

    //     switch (axis)
    //     {
    //         case 0:
    //             // yで積分
    //             break;
    //         case 1:
    //             // xで積分
    //             {
    //                 int i = 0;
    //                 for (auto& xvec : pdf)
    //                 {
    //                     T temp = 0.0f;
    //                     for (const auto& v : xvec)
    //                     {
    //                         temp += v;
    //                         // std::cout << "v = " << v << std::endl;
    //                         // std::cout << "temp = " << temp << std::endl;
    //                     }
    //                     // T temp = accumulate(xvec.begin(), xvec.end(), remove_reference<decltype(*std::begin(xvec))>::type{});
    //                     result[i++] = temp;
    //                 }
    //                 break;
    //             }
    //         default:
    //             break;
    //     }

    //     return result;
    // }

    // [0, 1]の範囲で乱数を生成
    float generate_random() const
    {
        random_device seed_gen;
        static std::mt19937 generator(seed_gen());
        static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
        return distribution(generator);
    }
    // Float generate_random() const
    // {
    //     random_device seed_gen;
    //     static std::mt19937 generator(seed_gen());
    //     static std::uniform_real_distribution<float> distribution(0.0f, 1.0f);
    //     return distribution(generator)*Float(1.0);
    // }

    sinc(const Properties &props) : Base(props)
    {
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
        
        m_divergence = props.texture<Texture>("divergence", 1.f);
        m_reflectance = props.texture<Texture>("reflectance", 1.f);
        // LUTFilename = props.set_string("LUTFilename", LUTFilename);
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

        // LUTFilename = props.texture<Texture>("LUTFilename", LUTFilename);
        m_a = props.texture<Texture>("a");
        // m_a = props.set_float("a", 1.f);
        // *m_a = 1.f;

        isTraversed = false;
        // LUTFilename = ""

        // ifstream ifs(LUTFilename.c_str());
        // LUTFilename = "/home/sugawara.ryo/BRDFEstimation/json/test.json";
        LUTFilename = "/home/sugawara.ryo/BRDFEstimation/pdf_sinc2, a = 10.00.json";
        std::cout << LUTFilename << std::endl;
        ifstream ifs(LUTFilename);
        if (ifs.good())
        {
            json m_json;
            ifs >> m_json;

            for (const auto& json : m_json)
            {
                // TODO floatで読んで、Tで入れるワンステップが必要かも
                vector<float> d_float = json.get<vector<float>>();
                // vector<T> d_Float(d_float.size());
                // int i = 0;
                // for (const auto& v : d_float)
                // {
                //     // d_Float[i++] = Float(1.0)*v;
                //     // d_float[i++] = v;
                // }
                // m_data.push_back(d_Float);
                m_data.push_back(d_float);
            }
            std:cout << "ファイルを読んだ" << LUTFilename << std::endl;
        }
        else
        {
            std::cout << "ファイルの読み込みに失敗しました" << std::endl;
        }

        // 周辺化を行う
        vector<float> marginalized_y = marginalize(m_data, 1);
        // vector<Float> marginalized_y = marginalize(m_data, 1);
        // yの累積分布関数を求める
        cum_theta = cumulative_sum(marginalized_y);
        // xごとの累積分布関数を求める
        
        std::cout << "json is read" << std::endl;
        std::cout << "test : " << m_data[1][5] << std::endl;
    }

    Vector3f rotate (const Vector3f vec, const Vector3f normal, const Vector3f to) const
    {
        Float cos_theta_i = Frame3f::cos_theta(to);
        Float sin_theta_i = Frame3f::sin_theta(to);
        Float cos_phi_i = Frame3f::cos_phi(to);
        Float sin_phi_i = Frame3f::sin_phi(to);

        Vector3f axis = dr::cross(normal, to);

        // // x軸周りの回転
        // Float temp_y = sin_theta_i*vec.y();
        // Float temp_z = cos_theta_i*vec.z();

        // // z軸周りの回転
        // Float result_x = cos_phi_i*vec.x();
        // Float result_y = sin_phi_i*temp_y;
        // Float result_z = temp_z;
        // ロドリゲスの回転公式
        Vector3f result = cos_theta_i*vec + dr::dot(axis, vec)*(1-cos_theta_i)*axis + dr::cross(axis, vec)*sin_theta_i;
        result = dr::normalize(result);
        // cout << result << endl;
        // return Vector3f(result_x, result_y, result_z);
        return result;
    }

    vector<vector<float>> make_sinc(vector<vector<float>>& vecs, int M, int N, float a, int expo) const
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
        v = normalize(vecs);
        for (int i = 0; i < N; i++) { 
            for (int j = 0; j < M; j++)
            {
                cout << "parts[" << i << "][" << j << "] = " << v[i][j] << endl;
            }
        }
        return v;
    }

    // vector<vector<Float>> make_sinc(vector<vector<Float>>& vecs, int M, int N, Float a_Float, int expo) const
    // {
    //     float a = a_Float[0];
    //     // vector<vector<float>> v;
    //     for (int i = 0; i < N; i++) { 
    //         vector<Float> parts = vector<Float>(M); // 0で初期化
    //         for (int j = 0; j < M; j++)
    //         {
    //             Float arg = a * static_cast<Float>(dr::sqrt( std::pow(j - M/2, 2) + std::pow(i - N/2, 2) ));
    //             parts[j] = arg == 0 ? 1.0f : dr::pow(dr::sin(arg), expo) / dr::pow(arg, expo); // argが0の時はゼロ除算を防ぐため
    //             cout << "parts["<< j << "] = " << parts[j] << endl;
    //         }
    //         vecs.push_back(parts);
    //         // v[i] = parts;
    //     }
    //     //正規化
    //     vector<vector<Float>> v;
    //     v = normalize(vecs);
    //     for (int i = 0; i < N; i++) { 
    //         for (int j = 0; j < M; j++)
    //         {
    //             cout << "parts[" << i << "][" << j << "] = " << v[i][j] << endl;
    //         }
    //     }
    //     return v;
    // }


    // template<class ForwardIterator>
    // ForwardIterator lower_bound_for_Float(ForwardIterator first, ForwardIterator last, const Float value) const
    // {
    //     using diff = typename std::iterator_traits<ForwardIterator>::difference_type;
    //     for (diff len = std::distance(first, last); len != 0; ) {
    //         diff half = len / 2;
    //         ForwardIterator mid = first;
    //         std::advance(mid, half);
    //         // if (*mid > value) {
    //         //     len -= half + 1;
    //         //     first = ++mid;
    //         // } else {
    //         //     len = half;
    //         // }
    //         // dr::if_stmt (
    //         Int lentemp = dr::select (*mid > value, len - half + 1, half);
    //         ForwardIterator firsttemp = dr::select (*mid > value, ++mid, first);
    //         len = lentemp;
    //         first = firsttemp;
    //     }
    //     return first;
    // }

    // template<class ForwardIterator>
    // ForwardIterator lower_bound_for_Float(ForwardIterator first, ForwardIterator last, const Float value) const
    // {
    //     for (Int8 len = std::distance(first, last)*Int8(1); len != 0; ) {
    //         Int8 half = len / 2;
    //         ForwardIterator mid = first;
    //         std::advance(mid, half);
    //         // if (*mid > value) {
    //         //     len -= half + 1;
    //         //     first = ++mid;
    //         // } else {
    //         //     len = half;
    //         // }
    //         // dr::if_stmt (
    //         len = dr::select (*mid > value, len - half + 1, half);
    //         // first = dr::select (*mid > value, ++mid, first);
    //     }
    //     return first;
    // }


    void traverse(TraversalCallback *callback) override
    {
        // std::string test_str = "a = "+std::to_string(0.1)+".json";
        // std::string test_str = "a = "+to_string(m_a.get())+".json";
        // std::string test_str = to_string();
        LUTFilename = make_LUTPath();
        std::cout << "LUTFilename = " << LUTFilename << std::endl;

        callback->put_object("divergence", m_divergence.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("reflectance", m_reflectance.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        // callback->put_parameter("a", *m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        callback->put_object("a", m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        // callback->put_parameter("a", *m_a.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
        // TODO LUTFilenameを加える？
        // callback->put_parameter("LUTFilename", LUTFilename, ParamFlags::NonDifferentiable);

        // std::cout << "Traverse is executing. a = " << m_a << std::endl;
        isTraversed = true;

        // mi.traverse(scene)が呼び出されたときに新しいaのLUTを作る
        if (isTraversed)
        {
            // a = m_a.get()[0];
            // m_data = make_sinc(m_data, 360, 1440, a, 2);
            // m_data = make_sinc(m_data, 360, 1440, m_a.get()[0], 2);
            // vector<vector<float>> data;
            // m_data = make_sinc(m_data, 360, 1440, 0.1f, 2);
            // m_data = data;

            // std::ostringstream out;
            // int precision = 2;
            // out << std::fixed << std::setprecision(precision) << a;
            // std::string a_str = out.str();
            // std::string filename = "pdf_sinc2, a = "+a_str+".json";

            // json m_json(m_data);

            // ofstream writing_file;
            // writing_file.open(filename, ios::out);
            // writing_file << m_json.dump(4) << endl;
            // writing_file.close();

            // isTraversed = false;


            // traverseの時にLoouup Tableをファイルから読み込む
            // ifstream ifs(LUTFilename.c_str());
            // LUTFilename = "/home/sugawara.ryo/BRDFEstimation/json/pdf_unitDist.json";
            // LUTFilename = "/home/sugawara.ryo/BRDFEstimation/json/pdf_sinc2, a = UniformSpectrum[value=0.100000].json";

            // もともとあるlUTを一旦削除
            for (auto& l : m_data)
            {
                l.clear();
            }
            m_data.clear();
            // 新たなLUTを読み込む
            ifstream ifs(LUTFilename);
            // ifstream ifs("pdf_sinc2, a = UniformSpectrum[value=0.100000].json");
            // ifstream ifs("test.json");
            if (ifs)
            {
                json m_json;
                ifs >> m_json;

                int j = 0;
                for (const auto& json : m_json)
                {
                    // std::cout << "Reading j = " << j << std::endl;
                    // TODO floatで読んで、Tで入れるワンステップが必要かも
                    vector<float> d_float = json.get<vector<float>>();
                    // vector<T> d_Float(d_float.size());
                    // int i = 0;
                    // for (const auto& v : d_float)
                    // {
                    //     // d_Float[i++] = Float(1.0)*v;
                    //     d_float[i++] = v;
                    // }
                    // m_data.push_back(d_Float);
                    m_data.push_back(d_float);
                    // m_data[j++] = d_float;
                }
                // std::cout << "File is read. The file path is " << LUTFilename << std::endl;
            }
            else
            {
                std::cout << "ファイルの読み込みに失敗しました" << std::endl;
            }

            // 周辺化を行う
            vector<float> marginalized_y = marginalize(m_data, 1);
            // vector<Float> marginalized_y = marginalize(m_data, 1);
            // yの累積分布関数を求める
            cum_theta = cumulative_sum(marginalized_y);
            // xごとの累積分布関数を求める
            
            std::cout << "json is read through traverse. " << LUTFilename << std::endl;
            std::cout << "test : " << m_data[1][5] << std::endl;
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

        // ----- 出射光のサンプリング -----
        // 最初にthetaのサンプリング
        auto itr_v_theta = lower_bound(cum_theta.begin(), cum_theta.end(), generate_random());
        // auto itr_u_phi = lower_bound_for_Float(cum_theta.begin(), cum_theta.end(), generate_random());
        // cout << *itr_u_phi << endl;
        auto index_theta = distance(cum_theta.begin(), itr_v_theta);
        // cout << index_phi << endl;
        // つづいてthetaのサンプリング　書き換え前
        // vector<float> pdf_theta = m_data[index_phi]; // thetaの分布を取得
        // vector<float> pdf_thetanorm = normalize(pdf_theta); // 正規化
        // vector<float> cum_x = cumulative_sum(pdf_thetanorm);
        // 書き換え後
        auto pdf_phi = m_data[index_theta]; // phiの分布を取得
        auto pdf_phinorm = normalize(pdf_phi); // 正規化
        auto cum_phi = cumulative_sum(pdf_phinorm);
        // auto itr_v_theta = lower_bound_for_Float(cum_x.begin(), cum_x.end(), generate_random());
        auto itr_u_phi = lower_bound(cum_phi.begin(), cum_phi.end(), generate_random());
        auto index_phi = distance(cum_phi.begin(), itr_u_phi);
        // cout << index_theta << endl;
        // cout << "cum_theta.size() = " << cum_theta.size() << endl;
        vector<Float> sampled_angles{index_phi*dr::TwoPi<Float>/cum_phi.size(), index_theta*dr::Pi<Float>/(2.0f*cum_theta.size())}; // サンプリングしたphi, thetaのペア
        // sampled_angles[0] -= dr::Pi<Float>;
        sampled_angles[0] = dr::TwoPi<Float> * generate_random(); // 角度の差分にする
        // sampled_angles[0] = 0.f;
        // sampled_angles[0] = dr::Pi<Float>/100.0f * generate_random(); // 角度の差分にする
        sampled_angles[1] -= dr::Pi<Float>/4.0f; // 角度の差分にする
        // sampled_angles[1] = 0.f; // 角度の差分にする
        // sampled_angles[1] = dr::Pi<Float>/2.0f * generate_random(); // 角度の差分にする

        // cout << sampled_angles << endl;
        Float del_phi = sampled_angles[0], del_theta = sampled_angles[1];        
        Float cos_del_phi = dr::cos(del_phi), sin_del_phi = dr::sin(del_phi);
        Float cos_del_theta = dr::cos(del_theta), sin_del_theta = dr::sin(del_theta);
        Vector3f delvec = Vector3f(sin_del_theta*cos_del_phi, sin_del_theta*sin_del_phi, cos_del_phi);
        Vector3f norm = Vector3f(0.f, 0.f, 1.f);
        Vector3f wo = rotate(delvec, norm, si.wi);
        bs.wo = wo;

        // Float cos_theta_o = cos_theta_i*cos_del_theta - sin_theta_i*sin_del_theta;
        // Float sin_theta_o = sin_theta_i*cos_del_theta + cos_theta_i*sin_del_theta;
        // Float cos_phi_o = cos_phi_i*cos_del_phi - sin_phi_i*sin_del_phi;
        // Float sin_phi_o = sin_phi_i*cos_del_phi + cos_phi_i*sin_del_phi;
        // Float wo_x = sin_theta_o*cos_phi_o;
        // Float wo_y = sin_theta_o*sin_phi_o;
        // Float wo_z = cos_theta_o;
        // Normal3f h = Vector3f(wo_x, wo_y, wo_z);

        // bs.wo = warp::square_to_cosine_hemisphere(sample2); // 出射光のサンプリング
        // Vector3f ideal_wo = reflect(si.wi);
        // bs.wo = Vector3f(wo_x, wo_y, wo_z) + ideal_wo;
        // bs.wo = Vector3f(wo_x, wo_y, wo_z);

        bs.sampled_component = 0;
        bs.sampled_type =+ BSDFFlags::DiffuseReflection;
        // bs.wo = 
        bs.eta = 1.f;
        bs.pdf = 1.f;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);

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

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active) * m_data[index_phi][index_theta];

        // return depolarizer<Spectrum>(value) & active;

        return 0.f;
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
        
        Float pdf = m_data[index_phi][index_theta];

        // return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
        return 0.f;
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

            active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

            UnpolarizedSpectrum value = m_reflectance->eval(si, active) * m_data[index_phi][index_theta];

            Float pdf = m_data[index_phi][index_theta];

            // return {depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f)};
            return {0.f, 0.f};
        }

    // std::string to_string() const override {
    //     std::ostringstream oss;
    //     oss << "sinc[" << std::endl
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
        // oss << "/home/sugawara.ryo/BRDFEstimation/json/pdf_sinc2, a = " << m_a.get() << ".json" << std::endl;
        std::ostringstream out;
        int precision = 2;
        out << std::fixed << std::setprecision(precision) << m_a.get()->max();
        std::string a_str = out.str();
        // oss << "/home/sugawara.ryo/BRDFEstimation/json/pdf_sinc2, a = " << a_str << ".json";
        oss << LUTDir << "pdf_sinc2, a = " << a_str << ".json";
        return oss.str();
    }


    MI_DECLARE_CLASS()

// フィールド
private:
    ref<Texture> m_divergence;
    ref<Texture> m_reflectance;
    // ref<float> m_a;
    // ref<Float> m_a;
    ref<Texture> m_a;
    vector<vector<float>> m_data;
    // vector<vector<Float>> m_data;
    vector<float> cum_theta;
    // vector<Float> cum_theta;
    std::string LUTFilename;
    std::string LUTDir;
    // String LUTFilename;
    // ref<Texture> LUTFilename;
    int index_phi, index_theta;

    bool isTraversed;
};

MI_IMPLEMENT_CLASS_VARIANT(sinc, BSDF)
MI_EXPORT_PLUGIN(sinc, "sinc")
NAMESPACE_END(mitsuba)