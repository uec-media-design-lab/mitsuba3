# Mitsuba Renderer 3 - AIRR
[Mitsuba Renderer 3]()のAIRRによる空中像シミュレーション用のブランチです。

以下のプラグインが新しく利用できます。

- 再帰反射材(roughretroreflector)
- ハーフミラー(halfmirror)

以下の既存のプラグインに変更を加えました。

- エリアライト(area)

基本的な使い方や他のプラグインの使い方などは[公式ドキュメント](https://mitsuba.readthedocs.io/en/stable/index.html)をご確認ください。

## 目次

- [Mitsuba Renderer 3 - AIRR](#mitsuba-renderer-3---airr)
  - [目次](#目次)
  - [チュートリアル](#チュートリアル)
  - [再帰性反射材(roughretroreflector)](#再帰性反射材roughretroreflector)
  - [ハーフミラー(halfmirror)](#ハーフミラーhalfmirror)
  - [エリアライト(area)](#エリアライトarea)
- [Mitsuba Renderer 3](#mitsuba-renderer-3)
  - [Introduction](#introduction)
  - [Main Features](#main-features)
  - [Tutorial videos, documentation](#tutorial-videos-documentation)
  - [Installation](#installation)
    - [Requirements](#requirements)
  - [Usage](#usage)
  - [About](#about)

## チュートリアル
[こちら](https://github.com/uec-media-design-lab/mitsuba3/blob/airr/airr-examples/render.ipynb)からAIRRのシーンを読み込み、それぞれの光学素子のパラメータを変えながら空中像をレンダリングする様子が確認できます。


## 再帰性反射材(roughretroreflector)

| Parameter | Type | Description |
| :---: | :---: | :---: |
| `int_ior` | `float` or `string` | プリズムの絶対屈折率を設定します。数値もしくは[マテリアル名](https://mitsuba.readthedocs.io/en/stable/src/generated/plugins_bsdfs.html#dielectric-ior-list:~:text=Instead%20of%20specifying%20numerical%20values%20for%20the%20indices%20of%20refraction%2C%20Mitsuba%203%20comes%20with%20a%20list%20of%20presets%20that%20can%20be%20specified%20with%20the%20material%20parameter%3A)で設定できます。(デフォルト：`bk7` / 1.5046) |
| `ext_ior` | `float` or `string` | 再帰反射材の外部の絶対屈折率を設定します。数値もしくは[マテリアル名](https://mitsuba.readthedocs.io/en/stable/src/generated/plugins_bsdfs.html#dielectric-ior-list:~:text=Instead%20of%20specifying%20numerical%20values%20for%20the%20indices%20of%20refraction%2C%20Mitsuba%203%20comes%20with%20a%20list%20of%20presets%20that%20can%20be%20specified%20with%20the%20material%20parameter%3A)から設定します。(デフォルト：`air` / 1.000277) |
| `base_material` | `string` | 裏面に蒸着された材質の屈折率を設定します。名前は[ドキュメント記載のマテリアル名](https://mitsuba.readthedocs.io/en/stable/src/generated/plugins_bsdfs.html#conductor-ior-list)から設定できます。(デフォルト：Al) |
| `base_eta`, `base_k` | `spectrum` or `texture` | 裏面に蒸着された材質の屈折率を設定します。`base_material`の代わりに複素屈折率で直接設定する場合は、この2つのパラメータを利用します。(デフォルト：`base_material`に従う) |
| `distribution` | `string` | 反射光の拡散に関連した微小表面の法線分布関数を設定します。`beckmann`と`ggx`が利用可能です。(デフォルト：`beckmann`)|
|`alpha`, `alpha_surface`, `alpha_internal`, `alpha_u_surface`, `alpha_v_surface`, `alpha_u_internal`, `alpha_v_internal` | `texture` or `float` | 各面の粗さを設定します。(デフォルト：`alpha`=0.1)<br> <hr> <br>`alpha`: すべての面の粗さを設定します。<br>`alpha_surface`: 表面のみの粗さを設定します。<br>`alpha_internal`: 裏面のみの粗さを設定します。<br>`alpha_u_surface`: 表面のUV座標系におけるU方向の粗さを設定します。<br>`alpha_v_surface`: 表面のUV座標系におけるV方向の粗さを設定します。<br>`alpha_u_internal`: 裏面のUV座標系におけるU方向の粗さを設定します。<br>`alpha_v_internal`: 裏面のUV座標系におけるV方向の粗さを設定します。|
| `surface_reflectance` | `spectrum` or `texture` | 表面の反射の係数を設定します。再帰反射における表面の透過にも影響します。(デフォルト：1.0) |
| `internal_reflectance` | `spectrum` or `texture` | 裏面の反射の係数を設定します。(デフォルト：1.0) |

実際のAIRRシステムに使われるような、裏面が金属蒸着されてできているコーナーキューブ型の再帰反射材のプラグインです。

`alpha`を設定することで、再帰反射の方向を中心に広がる様子を表現することができます。

```XML
<bsdf type="roughretroreflector">
    <float name="alpha" value="0.001"/>
</bsdf>

<bsdf type="roughretroreflector" id="retroreflector_smoothSurface">
    <float name="alpha_surface" value="0"/>
    <float name="alpha_internal" value="0.0005"/>
    <float name="internal_reflectance" vaue="0.83"/>
</bsdf>
```

## ハーフミラー(halfmirror)

| Parameter | Type | Description |
| :---: | :---: | :---: |
| `reflectance` | `float` | 反射率を設定します。(デフォルト：0.5) |
| `transmittance` | `float` | 透過率を設定します。(デフォルト：0.5) |

入射光を一定の確率で反射、一定の確率で透過する、ハーフミラーを模したプラグインです。

`reflectance`と`transmittance`の合計を1未満にすることで、ハーフミラーによる光の吸収も再現することができます。

```XML
<bsdf type="halfmirror" id="halfMirror_50R50T"/>

<bsdf type="halfmirror" id="halfMirror_55R30T">
    <float name="reflectance" value="0.55"/>
    <float name="transmittance" value="0.30"/>
</bsdf>
```


## エリアライト(area)

| Parameter | Type | Description |
| :---: | :---: | :---: |
| `radiance` | `spectrum` or `texture` | 放射輝度(単位面積・単位立体角あたりの放射束)を設定します。 |
| `coefficient` | `float` | `radiance`に一律に乗する係数を設定します。(デフォルト：1.0) |

Mitsuba Renderer 3に標準搭載されている[Area Light](https://mitsuba.readthedocs.io/en/stable/src/generated/plugins_emitters.html#area-light-area)に係数`coefficient`を加えたものです。

エリアライトの光は等方的に広がるため、どの方向から見ても見かけの明るさが等しくなります。

`radiance`に画像をテクスチャとして設定したうえで、`coefficient`のパラメータを調節することでディスプレイ自体の明るさをシーンファイルから設定できます。

```XML
<emitter type="area" id="display">
    <texture name="radiance" type="bitmap">
        <string name="filename" value="myDisplayedImage.png"/>
    </texture>
</emitter>

<emitter type="area" id="display_highBrightness">
    <texture name="radiance" type="bitmap">
        <string name="filename" value="myDisplayedImage.png"/>
    </texture>
    <float name="coefficient" value="2.0"/>
</emitter>
```


---

<!-- <img src="https://github.com/mitsuba-renderer/mitsuba3/raw/master/docs/images/logo_plain.png" width="120" height="120" alt="Mitsuba logo"> -->

<img src="https://raw.githubusercontent.com/mitsuba-renderer/mitsuba-data/master/docs/images/banners/banner_01.jpg"
alt="Mitsuba banner">

# Mitsuba Renderer 3

| Documentation  | Tutorial videos  | Linux             | MacOS             | Windows           |       PyPI        |
|      :---:     |      :---:       |       :---:       |       :---:       |       :---:       |       :---:       |
| [![docs][1]][2]| [![vids][9]][10] | [![rgl-ci][3]][4] | [![rgl-ci][5]][6] | [![rgl-ci][7]][8] | [![pypi][11]][12] |

[1]: https://readthedocs.org/projects/mitsuba/badge/?version=stable
[2]: https://mitsuba.readthedocs.io/en/stable/
[3]: https://rgl-ci.epfl.ch/app/rest/builds/buildType(id:Mitsuba3_LinuxAmd64Clang10)/statusIcon.svg
[4]: https://rgl-ci.epfl.ch/viewType.html?buildTypeId=Mitsuba3_LinuxAmd64Clang10&guest=1
[5]: https://rgl-ci.epfl.ch/app/rest/builds/buildType(id:Mitsuba3_LinuxAmd64gcc9)/statusIcon.svg
[6]: https://rgl-ci.epfl.ch/viewType.html?buildTypeId=Mitsuba3_LinuxAmd64gcc9&guest=1
[7]: https://rgl-ci.epfl.ch/app/rest/builds/buildType(id:Mitsuba3_WindowsAmd64msvc2020)/statusIcon.svg
[8]: https://rgl-ci.epfl.ch/viewType.html?buildTypeId=Mitsuba3_WindowsAmd64msvc2020&guest=1
[9]: https://img.shields.io/badge/YouTube-View-green?style=plastic&logo=youtube
[10]: https://www.youtube.com/watch?v=9Ja9buZx0Cs&list=PLI9y-85z_Po6da-pyTNGTns2n4fhpbLe5&index=1
[11]: https://img.shields.io/pypi/v/mitsuba.svg?color=green
[12]: https://pypi.org/pypi/mitsuba

## Introduction

Mitsuba 3 is a research-oriented rendering system for forward and inverse light
transport simulation developed at [EPFL](https://www.epfl.ch) in Switzerland.
It consists of a core library and a set of plugins that implement functionality
ranging from materials and light sources to complete rendering algorithms.

Mitsuba 3 is *retargetable*: this means that the underlying implementations and
data structures can transform to accomplish various different tasks. For
example, the same code can simulate both scalar (classic one-ray-at-a-time) RGB transport
or differential spectral transport on the GPU. This all builds on
[Dr.Jit](https://github.com/mitsuba-renderer/drjit), a specialized *just-in-time*
(JIT) compiler developed specifically for this project.

## Main Features

- **Cross-platform**: Mitsuba 3 has been tested on Linux (``x86_64``), macOS
  (``aarch64``, ``x86_64``), and Windows (``x86_64``).

- **High performance**: The underlying Dr.Jit compiler fuses rendering code
  into kernels that achieve state-of-the-art performance using
  an LLVM backend targeting the CPU and a CUDA/OptiX backend
  targeting NVIDIA GPUs with ray tracing hardware acceleration.

- **Python first**: Mitsuba 3 is deeply integrated with Python. Materials,
  textures, and even full rendering algorithms can be developed in Python,
  which the system JIT-compiles (and optionally differentiates) on the fly.
  This enables the experimentation needed for research in computer graphics and
  other disciplines.

- **Differentiation**: Mitsuba 3 is a differentiable renderer, meaning that it
  can compute derivatives of the entire simulation with respect to input
  parameters such as camera pose, geometry, BSDFs, textures, and volumes. It
  implements recent differentiable rendering algorithms developed at EPFL.

- **Spectral & Polarization**: Mitsuba 3 can be used as a monochromatic
  renderer, RGB-based renderer, or spectral renderer. Each variant can
  optionally account for the effects of polarization if desired.

## Tutorial videos, documentation

We've recorded several [YouTube videos][10] that provide a gentle introduction
Mitsuba 3 and Dr.Jit. Beyond this you can find complete Juypter notebooks
covering a variety of applications, how-to guides, and reference documentation
on [readthedocs][2].

## Installation

We provide pre-compiled binary wheels via PyPI. Installing Mitsuba this way is as simple as running

```bash
pip install mitsuba
```

on the command line. The Python package includes four variants by default:

- ``scalar_spectral``
- ``scalar_rgb``
- ``llvm_ad_rgb``
- ``cuda_ad_rgb``

The first two perform classic one-ray-at-a-time simulation using either a RGB
or spectral color representation, while the latter two can be used for inverse
rendering on the CPU or GPU. To access additional variants, you will need to
compile a custom version of Dr.Jit using CMake. Please see the
[documentation](https://mitsuba.readthedocs.io/en/latest/src/developer_guide/compiling.html)
for details on this.

### Requirements

- `Python >= 3.8`
- (optional) For computation on the GPU: `Nvidia driver >= 495.89`
- (optional) For vectorized / parallel computation on the CPU: `LLVM >= 11.1`

## Usage

Here is a simple "Hello World" example that shows how simple it is to render a
scene using Mitsuba 3 from Python:

```python
# Import the library using the alias "mi"
import mitsuba as mi
# Set the variant of the renderer
mi.set_variant('scalar_rgb')
# Load a scene
scene = mi.load_dict(mi.cornell_box())
# Render the scene
img = mi.render(scene)
# Write the rendered image to an EXR file
mi.Bitmap(img).write('cbox.exr')
```

Tutorials and example notebooks covering a variety of applications can be found
in the [documentation][2].

## About

This project was created by [Wenzel Jakob](https://rgl.epfl.ch/people/wjakob).
Significant features and/or improvements to the code were contributed by
[Sébastien Speierer](https://speierers.github.io/),
[Nicolas Roussel](https://github.com/njroussel),
[Merlin Nimier-David](https://merlin.nimierdavid.fr/),
[Delio Vicini](https://dvicini.github.io/),
[Tizian Zeltner](https://tizianzeltner.com/),
[Baptiste Nicolet](https://bnicolet.com/),
[Miguel Crespo](https://mcrespo.me/),
[Vincent Leroy](https://github.com/leroyvn), and
[Ziyi Zhang](https://github.com/ziyi-zhang).

When using Mitsuba 3 in academic projects, please cite:

```bibtex
@software{Mitsuba3,
    title = {Mitsuba 3 renderer},
    author = {Wenzel Jakob and Sébastien Speierer and Nicolas Roussel and Merlin Nimier-David and Delio Vicini and Tizian Zeltner and Baptiste Nicolet and Miguel Crespo and Vincent Leroy and Ziyi Zhang},
    note = {https://mitsuba-renderer.org},
    version = {3.1.1},
    year = 2022
}
```
