#!/usr/bin/python3

import mitsuba as mi
import xml.etree.ElementTree as ET
import math
import numpy as np


#spp = 4096
spp = 1

img_size = (1000, 1000)

float_dist = 50
cam_midair_dist = 200
display_z_offset = -25
lightsrc_scale = 30

elevation = 45
azimuth = 0

def setting(dist, azimuth, elevation):

    azimuth_rad = azimuth * math.pi / 180.0
    elevation_rad = elevation * math.pi / 180.0

    unit_pos = np.array([
        -math.sin(azimuth_rad) * math.cos(elevation_rad),
        math.sin(elevation_rad),
        -math.cos(azimuth_rad) * math.cos(elevation_rad)
    ])

    display_offset = np.array([
        display_z_offset * math.sin(azimuth_rad),
        0,
        display_z_offset * math.cos(azimuth_rad)
    ])

    display_pos = unit_pos * dist + display_offset
    cam_pos = unit_pos * cam_midair_dist + display_pos

    return display_pos, cam_pos


def save_xml(tree, display_pos, azimuth, elevation, fileName):

    root = tree.getroot()
    
    updates = {
        'lightsrc_scale': str(lightsrc_scale),
        'lightsrc_pos_x': str(display_pos[0]),
        'lightsrc_pos_y': str(-display_pos[1]),
        'lightsrc_pos_z': str(display_pos[2]),
        'lightsrc_rot_x': str(-elevation),
        'lightsrc_rot_y': str(azimuth),
    }

    for default in root.findall('default'):
        name = default.get('name')
        if name in updates:
            default.set('value', updates[name])

    tree.write(fileName, encoding='utf-8', xml_declaration=False)


def get_cam(cam_pos, display_pos):

    cam = mi.load_dict({
        # typeフィールド：作成するオブジェクトの名前を定義する
        'type': 'perspective',
        # その他フィールド：そのオブジェクトのプロパティを指定する
        'to_world': mi.ScalarTransform4f.look_at(
            target = [display_pos[0], display_pos[1], display_pos[2]],
            origin = [cam_pos[0], cam_pos[1], cam_pos[2]],
            up = [0, 1, 0]
        ),
        # フィルムインスタンスの辞書をネストしている
        'film': {
            'type': 'hdrfilm',
            'width': img_size[0], 'height':img_size[1],
        },
        # focal_lengthの代わりにfovを設定する（焦点距離50mm、フィルムサイズ：35mmフル）
        'fov_axis': 'x',
        'fov': 39.5978
    })

    return cam


tree = ET.parse('model.xml')

xml_filename = "{:04}_{:04}.xml".format(elevation, azimuth)
display_pos, cam_pos = setting(float_dist, azimuth, elevation)
save_xml(tree, display_pos, azimuth, elevation, xml_filename)

mi.set_variant('cuda_ad_rgb')

cam = get_cam(cam_pos, display_pos)
scene = mi.load_file(xml_filename)
image = mi.render(scene, sensor=cam, spp=spp)
mi.util.write_bitmap("result.png", image)
