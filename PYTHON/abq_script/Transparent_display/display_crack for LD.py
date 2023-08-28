#!/user/bin/python
# -*- coding: utf-8 -*-

#=======================利用显示组设置显示裂纹===================================
#liuxin66@live.com
#2020年7月14日
# V2.0
# 运行方式：
#   1、菜单栏plug-ins...
#===============================================================================
from odbAccess import *
from  abaqus import *
from abaqusConstants import *
import displayGroupOdbToolset
import sys
sys.path.append('G:\\CAE\\Python\\预制裂纹')  #填写当前"facelements_L.py 、facelements_w.py"目录
#from facelements_L import *
from facelements_w import *
###############集合定义###############
#region
   
def exmaple():
    elementLabels_1 = (

        )
    elementLabels_2 = (

        )
    elementLabels_3 = (

        )
    elementLabels_4 = (

        )
    return elementLabels_1,elementLabels_2,elementLabels_3,elementLabels_4

def create_surfaces(instName,surfName,elementLabels,odb1):
    a    = odb1.rootAssembly
    ( e_l1 ,e_l2 ,e_l3 ,e_l4 ) = elementLabels

    try:
        e_set1 = a.ElementSetFromElementLabels(name = 'e_set1' , elementLabels = ( ( instName, e_l1 ), ) )
        e_set2 = a.ElementSetFromElementLabels(name = 'e_set2' , elementLabels = ( ( instName, e_l2 ), ) )
        e_set3 = a.ElementSetFromElementLabels(name = 'e_set3' , elementLabels = ( ( instName, e_l3 ), ) )
        e_set4 = a.ElementSetFromElementLabels(name = 'e_set4' , elementLabels = ( ( instName, e_l4 ), ) )
    except:
        print "e_set* has been defined by others"
        e_set1 = a.elementSets['e_set1']
        e_set2 = a.elementSets['e_set2']
        e_set3 = a.elementSets['e_set3']
        e_set4 = a.elementSets['e_set4']
    else:
        print 'Create SETs: e_set1,  e_set2,   e_set3,    e_set4'
    
    elementSetSeq = ( (e_set1, FACE1), (e_set2, FACE2), (e_set3, FACE3) , (e_set4, FACE4) )

    try:
        a.MeshSurfaceFromElsets(name=surfName, elementSetSeq=elementSetSeq)               #创建surf集合
    except:
        print "'OUTSIDE_SURFS' has been defined by others"
    else:
        print 'OUTSIDE_SURFS has been created'
    
def create (odb1):
    # input "o" session.odbs
    # create surface : OUTSIDE_SURFS

    name = odb1.name
    instName    = 'PART-1-1'
    surfName    = 'OUTSIDE_SURFS'
    
    # 需要修改ODB的模型名称，在如下的if条件判断语句中！！！
    if   'Job_Test4_120MPa_5' in name :       
        create_surfaces(    instName,    surfName,    try_confining(),   odb1   )
    elif 'LD15_V0607' in name :
        create_surfaces(    instName,    surfName,    LD15(),   odb1   )
    elif 'LD10_V0607' in name :
        create_surfaces(    instName,    surfName,    LD10_AL(),   odb1   )
    elif 'LD25CS0_v0607' in name :
        create_surfaces(    instName,    surfName,    LD25_AL(),   odb1   )
    elif 'LD20_V0607' in name :
        create_surfaces(    instName,    surfName,    LD20(),   odb1   )
    elif 'W1' in name :
        create_surfaces(    instName,    surfName,    Fissure30_W1(),   odb1   )
    elif 'W3' in name :
        create_surfaces(    instName,    surfName,    Fissure30_W3(),   odb1   )            
    else:
        print '模型没有被添加，请添加该模型的单元集合'
#endregion
#####################################

#================================================输入================================================                                               #新建窗口的名字
solid_elset_name  ='PART-1-1.C3D10M'
coh_elset_name    ='COH3D6'
input_parameter   =(solid_elset_name,coh_elset_name)
#path              ='G:/CAE/BDS-25-wedge-v4.odb'
#================================================裂纹显示函数=========================================
def display_crack3(i,s_odb,s_view,input_parameter):
    '''
    以显示coh单元的方式显示裂纹，coh单元只显示上下表面
    '''
    instName    = 'PART-1-1'                                        #岩石部件名称
    setName     = 'COH3D6' 
    solid_elset_name    = input_parameter[0]
    coh_elset_name      = input_parameter[1]
    odb_num = s_odb[i]

    view_name   ='Viewport: ' +str ( s_view[i] )
    dgo         = displayGroupOdbToolset  

    os1 = session.odbs.values()[odb_num]                            #指定ODB
    myAssembly = os1.rootAssembly
    inst    = myAssembly.instances[instName]
    coh_set = inst.elementSets[setName]                             #指定set
    step1   = os1.steps.values()[0] 
    #region 获取最后一个分析步的裂纹情况，确定裂纹的统计区域(只对会断裂的coh单元进行统计提高效率)
    #----------------------------------------------------
    #创建裂纹区域Set
    c_name = 'crack_elements_finalframe'

    try:
        crack_coh_set = inst.elementSets[c_name]
    except KeyError:
        try:
            final_frame_sdv = step1.frames[-1].fieldOutputs['SDV1'].getSubset(region = coh_set).values
        except :
            final_frame_sdv = step1.frames[-1].fieldOutputs['SDV_STATUS'].getSubset(region = coh_set).values

        final_sdv_list  = []
        final_ele_lab_list = []
        for i in final_frame_sdv:                   #将sdv值按顺序存储到列表中，便于处理
            final_sdv_list.append(i.data)
            final_ele_lab_list.append(i.elementLabel)

        f_sdv_ip_1  = final_sdv_list [0::3]         #将sdv_element每3个元素中的第一个取出构成新数组（按积分点分成3个列表）
        f_sdv_ip_2  = final_sdv_list [1::3]
        f_sdv_ip_3  = final_sdv_list [2::3]
        f_ele_lab   = final_ele_lab_list[0::3]
        max_crack   = len(f_sdv_ip_1)


        lab_crack_ele     = []                      #统计满足条件的  单元编号
        for i in range( 0,max_crack ):              #遍历所有coh单元
            status = f_sdv_ip_1[i] + f_sdv_ip_2[i] + f_sdv_ip_3[i]   

            if status == 0:                         #当有单元发生破坏时

                lab_crack_ele.append(f_ele_lab[i])  #统计满足条件的 evol列表的 单元编号

        #为满足条件的单元创建集合
        inst.ElementSetFromElementLabels(name = c_name, elementLabels=lab_crack_ele)
        crack_coh_set = inst.elementSets[c_name]
        print "'%s' has been creacted" %c_name
        print 'MAX %s cohsive is cracked ' %max_crack
    else:  
        print "'%s' has been defined by others. \nPlease check that the SET is correct,"%c_name\
        +"or reopen the ODBs to create it"
    #----------------------------------------------------
    #endregion

    cohsurfs_name    = 'COH_SURFS'
        
    session.Viewport(name=view_name, origin=(20, 20), width=160, height=100)            #创建新窗口
    session.viewports[view_name].setValues(displayedObject=os1)                         #在窗口中显示odb
    view1 = session.viewports[view_name].odbDisplay                                     #赋值，缩短命令
    view1.setStatusVariable(variableLabel='SDV1', outputPosition=INTEGRATION_POINT, useStatus=True, 
        statusMinimum=0.5)                                                              #设置仅显示已删除的coh单元


    #创建coh表面
    e_set = os1.rootAssembly.instances['PART-1-1'].elementSets[c_name]          
    elementSetSeq = ((e_set,FACE1),(e_set,FACE2))
    try:
        os1.rootAssembly.MeshSurfaceFromElsets(name=cohsurfs_name, elementSetSeq=elementSetSeq)#创建coh上下表面集合
    except:
        print "'%s' has been defined by others" %cohsurfs_name
    else:
        print 'Create surfaces:',cohsurfs_name

    #选择配色方案
    session.viewports[view_name].enableMultipleColors()
    session.viewports[view_name].setColor(initialColor='#BDBDBD')
    cmap=session.viewports[view_name].colorMappings['Element type']
    cmap.updateOverrides(overrides={'COH3D6':(True, '#333333', 'Default', 
    '#333333')})
    session.viewports[view_name].setColor(colorMapping=cmap)
    session.viewports[view_name].disableMultipleColors()


    #创建显示组（display group）                                                    

    leaf1 = dgo.LeafFromSurfaceSets(surfaceSets=(cohsurfs_name, ))                         #由单元set创建leaf
    session.DisplayGroup(name='DisplayGroup-crack', leaf=leaf1)                         #创建显示组
    dg_c  = session.displayGroups['DisplayGroup-crack']

    leaf2 = dgo.LeafFromElementSets(elementSets=(solid_elset_name)) 
    session.DisplayGroup(name='DisplayGroup-solid', leaf=leaf2)
    dg_s  = session.displayGroups['DisplayGroup-solid']

    view1.setValues(visibleDisplayGroups=(dg_c,dg_s))                                   #plot display group

    view1.display.setValues(plotState=(DEFORMED))

    view1.commonOptions.setValues(visibleEdges=FREE,edgeLineThickness=THIN,             #cohsive
                                translucencyFactor=0.05,translucency=OFF)

    view1.displayGroupInstances['DisplayGroup-crack'].setValues(lockOptions=ON)         #solid

    view1.commonOptions.setValues(visibleEdges=NONE,translucency=OFF)

    if len(s_odb)==1 :
        session.viewports[view_name].maximize()
    else:
        session.viewports[view_name].restore()
def display_crack2(i,s_odb,s_view,input_parameter):
    '''
    显示coh上下表面的方式显示裂纹；同时实体单元只显示模型表面
    最终可以实现透明视角观察裂纹
    需要建模时设定好表面集合（推荐） 或者 在odb中重建
    '''
    coh_elset_name      = input_parameter[1]
    odb_num = s_odb[i]                          #odb 索引号

    view_name         ='Viewport: ' +str ( s_view[i] )
    dgo   = displayGroupOdbToolset  

    os1 = session.odbs.values()[odb_num]                                                #指定ODB

    cohsurfs_name    = 'COH_SURFS'
    solidsurfs_name  = 'OUTSIDE_SURFS'

    #创建coh表面
    e_set = os1.rootAssembly.instances['PART-1-1'].elementSets[coh_elset_name]          
    elementSetSeq = ((e_set,FACE1),(e_set,FACE2))
    try:
        os1.rootAssembly.MeshSurfaceFromElsets(name=cohsurfs_name, elementSetSeq=elementSetSeq)#创建coh上下表面集合
    except:
        print "'%s' has been defined by others" %cohsurfs_name
    else:
        print 'Create surfaces:',cohsurfs_name

    #创建solid外表面
    create(os1)
        
    session.Viewport(name=view_name, origin=(20, 20), width=160, height=100)            #创建新窗口
    session.viewports[view_name].setValues(displayedObject=os1)                         #在窗口中显示odb
    view1 = session.viewports[view_name].odbDisplay                                     #赋值，缩短命令
    view1.setStatusVariable(variableLabel='SDV1', outputPosition=INTEGRATION_POINT, useStatus=True, 
        statusMinimum=0.5)                                                              #设置仅显示已删除的coh单元

    #创建显示组（display group）
                                                     
    #coh显示组
    leaf1 = dgo.LeafFromSurfaceSets(surfaceSets=(cohsurfs_name, ))                      #由单元set创建leaf
    session.DisplayGroup(name='DisplayGroup-crack', leaf=leaf1)                         #创建显示组
    dg_c  = session.displayGroups['DisplayGroup-crack']

    #solid显示组
    leaf2 = dgo.LeafFromSurfaceSets(surfaceSets=(solidsurfs_name, )) 
    session.DisplayGroup(name='DisplayGroup-solid', leaf=leaf2)
    dg_s  = session.displayGroups['DisplayGroup-solid']

    #刀具显示组
    leaf3 = dgo.LeafFromElementSets(elementSets=("CUTTER-1.SET-CELL", 
    "CUTTER-2.SET-CELL", ))
    session.DisplayGroup(name='DisplayGroup-cutter', leaf=leaf3)
    dg_ct  = session.displayGroups['DisplayGroup-cutter']

    view1.setValues( visibleDisplayGroups=(dg_c,dg_s,dg_ct) )                           #plot display group

    view1.display.setValues(plotState=(DEFORMED))
    
    view1.commonOptions.setValues(visibleEdges=FREE,edgeLineThickness=VERY_THIN,        #COH显示方式
                                edgeLineStyle=SOLID,translucencyFactor=0.25,translucency=ON)
    view1.displayGroupInstances['DisplayGroup-crack'].setValues(lockOptions=ON)

    view1.commonOptions.setValues(visibleEdges=NONE,translucency=ON,                    #solid显示方式
                                   translucencyFactor=0.3 ) 
    view1.displayGroupInstances['DisplayGroup-solid'].setValues(lockOptions=ON)

    view1.commonOptions.setValues(translucency=ON,edgeColorFillShade='#333333',         #cutter显示方式
    translucencyFactor=0.15,visibleEdges=FEATURE,renderStyle=SHADED)
    view1.displayGroupInstances['DisplayGroup-cutter'].setValues(lockOptions=ON)

    view1.commonOptions.setValues(translucency=OFF)                                     #关闭透明度避免显示卡慢    

    #选择配色方案
    session.viewports[view_name].enableMultipleColors()
    session.viewports[view_name].setColor(initialColor='#BDBDBD')
    cmap=session.viewports[view_name].colorMappings['Element type']
    cmap.updateOverrides(overrides={'COH3D6':(True, '#FF0000', 'Default',               #coh颜色
    '#FF0000')})
    cmap.updateOverrides(overrides={'C3D8':(True, '#008080', 'Default',                 #cutter颜色
    '#008080')})
    cmap.updateOverrides(overrides={'C3D10M':(True, '#F2F6D1', 'Default',               #岩石solid颜色
    '#F2F6D1')})
    session.viewports[view_name].setColor(colorMapping=cmap)
    session.viewports[view_name].disableMultipleColors()

    #调整光线
    session.viewports[view_name].lightOptions.lights[0].setValues(
    diffuseColor='#4d4d4d',specularColor='#878787',latitude=29,longitude=9)
    session.viewports[view_name].lightOptions.lights[1].setValues(
    diffuseColor='#424242',specularColor='#b5b5b5',latitude=33,longitude=-52)
    session.viewports[view_name].lightOptions.lights[2].setValues(
    diffuseColor='#2e2e2e',specularColor='#7a7a7a',latitude=-44,longitude=-4)
    session.viewports[view_name].lightOptions.setValues(ambientColor='#707070')


    if len(s_odb)==1 :
        session.viewports[view_name].maximize()
    else:
        session.viewports[view_name].restore()
def display_crack1(i,s_odb,s_view,input_parameter):
    '''
    以显示coh单元的方式显示裂纹
    '''
    solid_elset_name    = input_parameter[0]
    coh_elset_name      = 'PART-1-1.'+input_parameter[1]
    odb_num = s_odb[i]

    view_name   ='Viewport: ' +str ( s_view[i] )
    dgo         = displayGroupOdbToolset  

    os1 = session.odbs.values()[odb_num]                                                #指定CAE界面已打开的第一个odb
        
    session.Viewport(name=view_name, origin=(20, 20), width=160, height=100)            #创建新窗口
    session.viewports[view_name].setValues(displayedObject=os1)                         #在窗口中显示odb
    view1 = session.viewports[view_name].odbDisplay                                     #赋值，缩短命令

    try:
        view1.setStatusVariable(variableLabel='SDV1', outputPosition=INTEGRATION_POINT, useStatus=True, 
            statusMinimum=0.5)                                                              #设置仅显示已删除的coh单元
    except:
        view1.setStatusVariable(variableLabel='SDV_STATUS', outputPosition=INTEGRATION_POINT, useStatus=True, 
            statusMinimum=0.5)  

    #选择配色方案
    session.viewports[view_name].enableMultipleColors()
    session.viewports[view_name].setColor(initialColor='#BDBDBD')
    cmap=session.viewports[view_name].colorMappings['Set']
    cmap.updateOverrides(overrides={
        'PART-1-1.C3D10M':(True, '#F5F5DC', 'Default','#F5F5DC'),
        'PART-PLATE-1.SET-ALL':(True, '#448A73', 'Default', '#448A73'),
        'PART-INDENTER-1.SET-ALL':(True, '#9370DB', 'Default', '#9370DB'),
        'PART-1-1.COH3D6':(True, '#EE0000', 'Default', '#EE0000'),
        })
    session.viewports[view_name].setColor(colorMapping=cmap)
    session.viewports[view_name].disableMultipleColors()


    #创建显示组（display group）                                                    

    leaf1 = dgo.LeafFromElementSets(elementSets= (coh_elset_name))                      #由单元set创建leaf
    session.DisplayGroup(name='DisplayGroup-crack', leaf=leaf1)                         #创建显示组
    dg_c  = session.displayGroups['DisplayGroup-crack']

    leaf2 = dgo.LeafFromElementSets(elementSets=(solid_elset_name)) 
    session.DisplayGroup(name='DisplayGroup-solid', leaf=leaf2)
    dg_s  = session.displayGroups['DisplayGroup-solid']

    leaf3 = dgo.LeafFromElementSets(elementSets=('PART-INDENTER-1.SET-ALL',"PART-PLATE-1.SET-ALL", )) 
    session.DisplayGroup(name='DisplayGroup-indenter', leaf=leaf3)
    dg_i  = session.displayGroups['DisplayGroup-indenter']

    view1.setValues(visibleDisplayGroups=(dg_c,dg_s,dg_i))                                   #plot display group

    view1.display.setValues(plotState=(DEFORMED))

    view1.commonOptions.setValues(visibleEdges=FREE,edgeLineThickness=THIN,             #cohsive
                                translucencyFactor=0.05,translucency=ON)

    view1.displayGroupInstances['DisplayGroup-crack'].setValues(lockOptions=ON)         #锁定

    view1.commonOptions.setValues(visibleEdges=NONE,translucency=OFF)                   #solid

    if len(s_odb)==1 :
        session.viewports[view_name].maximize()
    else:
        session.viewports[view_name].restore()
#================================================选择函数=============================================
#获取要处理的模型范围

def simple (name):
    '''
    从odbname中去掉路径，得到文件名
    '''
    if '/' in name :                                       
        simple_name = name.split('/')[-1].split('.')[0]
    else:
        simple_name = name
    return  simple_name

def Select_object(name,objects):
    '''
    选择要处理的对象，和处理顺序
    返回object、Viewport的索引值列表\n
    使用示例：
           name ='JOBs'
           obj = mdb.jobs.values()
           (joblist,viewlist) = Select_object(name,obj)
           如输入为空则全选所有object
    '''
    title   = 'Select the %s from below\n\n'  %name
    count   = 0
    #输出job列表
    for i in objects:
        i_name = simple(i.name)

        if '/' in i_name :                                       
            i_name = i_name.split('/')[-1].split('.')[0]

        title += str(count) + '---' + i_name +'\n' 
        count += 1

    title += 'Example:  [0,2,4]  0,2,4  will be slected\n\
                [1:3]    1,2,3   will be slected\n               Multiselect  has a higher priority'
    fields = (('Multiselect mode',''),('Beginninig of slice',''),('End of slice',''),('Viewport','Default'),('Method:1,2,3','3' ))

    list_odject,slice_a,slice_b,list_Viewport,method = getInputs ( fields=fields,label=title )


    #获得object排序列表
    if list_odject == '':                               #如果为空则按slice生成列表
        if slice_a == '' and slice_b == '':
            s_object = range(count)
        else:
            slice_a = int (slice_a)
            slice_b = int (slice_b) + 1

            s_object   = range(slice_a,slice_b,1)
    elif list_odject =='all':
        s_object = range(count)
    else:
        s_object = map(int,list_odject.split(','))      #object 从0开始编号

    num_s_object = range(len(s_object))

    #获得Viewport排序列表
    if list_Viewport =='Default':
        s_Viewport = [i+1 for i in num_s_object]        #Viewport 从1开始编号
    else:
        s_Viewport = map(int,list_Viewport.split(','))

    print'\n-----------------------------'
    print"Order of Data ："
    for i in num_s_object:
        print '\t',i+1,simple(objects[i].name)
    print'\n-----------------------------'

    #method 
    method = int(method) 

    return s_object,s_Viewport,method

#=====================================================================================================

#=====main=====
def main():
    obj  = session.odbs.values()
    name = 'Odb'
    (odb_list,view_list,method) = Select_object(name,obj)

    if method== 1 :
        for  i in range(len(odb_list)):
            display_crack1(i,odb_list,view_list,input_parameter)
    elif method== 2:
        for  i in range(len(odb_list)):
            display_crack2(i,odb_list,view_list,input_parameter)
    elif method== 3:
        for  i in range(len(odb_list)):
            display_crack3(i,odb_list,view_list,input_parameter)    
    else:
        print "please input '1' or '2'"
#==============

if __name__ =='__main__':
    print   '\n//////////作为主程序调用//////////'
    main()  #作为主程序调用时执行
else:
    print '\n\n---------------------%s被其他程序调用---------------------'%__file__.split("\\")[-1]