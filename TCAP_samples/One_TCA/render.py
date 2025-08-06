#!/usr/bin/python3

import mitsuba as mi

#mi.set_variant('scalar_rgb')
#mi.set_variant('llvm_ad_rgb')
mi.set_variant('cuda_ad_rgb')

scene = mi.load_file("./model.xml")
spp = 1024

params = mi.traverse(scene)
#print(params)

cam1 = mi.load_dict({
    # typeフィールド：作成するオブジェクトの名前を定義する
    'type': 'perspective',
    # その他フィールド：そのオブジェクトのプロパティを指定する
    'to_world': mi.ScalarTransform4f.look_at(
        target = [0, 0, 0],
        origin = [-10, 2, 2],
        up = [0, 1, 0]
    ),
    # フィルムインスタンスの辞書をネストしている
    'film': {
        'type': 'hdrfilm',
        'width': 500, 'height':500,
    }
})

image = mi.render(scene, sensor=cam1, spp=spp)

mi.util.write_bitmap("one_tca.png", image)
