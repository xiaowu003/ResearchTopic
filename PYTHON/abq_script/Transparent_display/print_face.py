#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from abaqus import *
from abaqusConstants import *

def part_name(name):
    # 本函数用于生成没跑一次原模型提取的surface的表头
    m_name = file(output_path, 'a')
    m_name.write('def ' + p_name + '():\n')

def print_setlable(set_name, set_list, line_n, output_path):
    # type: (四个面列表名称, 面列表, 切分后每行的个数，输出的py文件) -> object
    # 用于导出odb或者cae文件中某个set的单元列表，输出相应py文件
    # examaple: set_elements = myAssembly.instances[instName].elementSets['name'].elements
    # line_n 每一行写入的编号数目
    # output_path 输出路径


    # 列表长度情况
    in_s = len(set_list) // line_n
    re_s = len(set_list) % line_n
    line_nn = range(0, line_n)
    # 写入数据
    f = file(output_path, 'a')
    # f.write()
    # f.write('def '+ model_name + '():' '\n')
    f.write('\t' + set_name + ' =(\n')


    for i in range(0, in_s):  # 写入完整行
        data_line = '\t\t'
        for n in line_nn:
            data_line += '%d,' % (set_list[n + i * line_n])
        data_line += '\n'
        f.write(data_line)

    data_line = '\t\t'
    for m in range(0, re_s):  # 写入非完整行
        data_line += '%d,' % (set_list[in_s * line_n + m])
    data_line += '\n'
    f.write(data_line)

    f.write('\t\t)\n')  # 结束
    f.close()

def return_face():
    # 用于表末写出 return
    m_name = file(output_path, 'a')
    m_name.write('\treturn face1elements,face2elements,face3elements,face4elements\n')

# 在part模块下建立surface
p = mdb.models['try'].parts['PART-1']          # 不同模型需要修改名称！！！
p_name = 'try_confining'                                    # 导出函数名称

surf1 = p.allSurfaces['Surf-1']                          #Surf-1为需要透明显示的外表面
surf1_ele = surf1.elements
surf1_side = surf1.sides
# 将由四边形单元四个face建立的surface分别输出
face1elements = []
face2elements = []
face3elements = []
face4elements = []


for i in range(len(surf1_side)):
    if surf1_side[i] == FACE1:
        face1elements.append(surf1_ele[i].label)

    elif surf1_side[i] == FACE2:
        face2elements.append(surf1_ele[i].label)

    elif surf1_side[i] == FACE3:
        face3elements.append(surf1_ele[i].label)

    else:
        face4elements.append(surf1_ele[i].label)

# p_key = mdb.models.keys()
# # print p_key[0]
# p_name = p_key[0]
# for i in range(len(p_name)):
#     if p_name[i] == '-':
#         p_name = p_name[:i] + p_name[i+1:]
# # p_name.replace('-', '_')
# print p_name

face1 = 'face1elements'
face2 = 'face2elements'
face3 = 'face3elements'
face4 = 'face4elements'
num = 10
output_path = 'facelements_w.py'
part_name(p_name)


print_setlable(face1, face1elements, num, output_path)
print_setlable(face2, face2elements, num, output_path)
print_setlable(face3, face3elements, num, output_path)
print_setlable(face4, face4elements, num, output_path)

return_face()


# return face1elements,face2elements,face3elements,face4elements

print('Write success!'+ ',data in ' + output_path)

# print face4elements
# print len(face4elements)



