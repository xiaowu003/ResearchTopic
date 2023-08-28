#!/user/bin/python
# -*- coding: utf-8 -*-
#liuxin66@live.com
# 统计coh裂纹面积
# 2020年7月16日 
# 运行方式：
#   1、菜单栏File->Run Script...
#================================================================================
from abaqus import *
from abaqusConstants import *
from odbAccess import *
import time
#==============================================函数定义==============================================
def get_crack(odb_num,input_parameter):
    output_path = input_parameter[0] 
    instName    = input_parameter[1]                                       
    setName     = input_parameter[2] 

    odb = session.odbs.values()[odb_num]                            #选择已打开的odb中的第一个
    myAssembly = odb.rootAssembly
    inst    = myAssembly.instances[instName]
    coh_set = inst.elementSets[setName]                             #指定set
    step1   = odb.steps.values()[0]                                 #指定分析步
    print "\n==================started=================="
    odb_name = odb.name.split('/')[-1].split('.')[0]
    l_name  = odb_name.split('-')
    xy_name = l_name[1]+'-'+l_name[-1]                               #xydata 数据后缀
    print xy_name
    print 'Get data from',odb_name

    #region 获取最后一个分析步的裂纹情况，确定裂纹的统计区域(只对会断裂的coh单元进行统计提高效率)
    #----------------------------------------------------
    final_frame_sdv = step1.frames[-1].fieldOutputs['SDV1'].getSubset(region = coh_set).values
    final_sdv_list  = []
    final_ele_lab_list = []
    for i in final_frame_sdv:                   #将sdv值按顺序存储到列表中，便于处理
        final_sdv_list.append(i.data)
        #final_ele_lab_list.append(i.elementLabel)

    f_sdv_ip_1  = final_sdv_list [0::3]         #将sdv_element每3个元素中的第一个取出构成新数组（按积分点分成3个列表）
    f_sdv_ip_2  = final_sdv_list [1::3]
    f_sdv_ip_3  = final_sdv_list [2::3]
    #f_ele_lab   = final_ele_lab_list[0::3]
    max_crack   = len(f_sdv_ip_1)


    number_crack_ip   = []
    number_crack_ele  = []                      #统计满足条件的 evol列表的 序号
    #lab_crack_ele     = []                      #统计满足条件的 evol列表的 单元编号
    for i in range( 0,max_crack ):              #遍历所有coh单元
        status = f_sdv_ip_1[i] + f_sdv_ip_2[i] + f_sdv_ip_3[i]   

        if status < 3:                          #当有单元发生破坏时

            number_crack_ele.append(i)          #统计满足条件的 evol列表的 序号
            #lab_crack_ele.append(f_ele_lab[i])  #统计满足条件的 evol列表的 单元编号

            n=3*i             
            number_crack_ip.extend([n,n+1,n+2]) #统计满足条件的 sdv1列表的 序号
            
    print 'MAX %s cohsive is cracked ' %max_crack
    #----------------------------------------------------
    #endregion

    #region 循环求取不同时刻的体积
    #----------------------------------------------------
    totalTime = 0.0                                                 #记录分析步时间
    count=0                                                         #记录次数
    frame_time_list=[]
    crack_length_list=[]
    crack_data=[]
    evol_element = step1.frames[0].fieldOutputs['EVOL'].getSubset(region = coh_set).values
    evol_ele_list =[]
    # evol_ele_lab_list =[]                                         #输出对应的element（调试用）
    for i in number_crack_ele:                                      #将数据转换成list
        evol_ele_list.append(evol_element[i].data)
        # evol_ele_lab_list.append(i.elementLabel)

    #frame= step1.frames[2]
    #指定evol的帧（所有帧内都一样，所以只输一次提高效率）
    #for step in odb.steps.values():                                #遍历所有分析步
    for frame in step1.frames:                                      #遍历所有帧
            print step1.name,'frame=',count,'...'
            count += 1
            crack_length = 0                                        #每帧裂纹长度赋初值
            frameTime   = frame.frameValue + totalTime              #当前帧对应的总时间
            
            sdv_element   = frame.fieldOutputs['SDV1'].getSubset(region = coh_set).values #sdv场变量输出的是积分点上的值，一个coh3d6有3个积分点

            #将数据转换成list
            sdv_ele_list        = []
            # sdv_ele_lab_list    = []                              #输出对应的element（调试用）
            # sdv_ele_ipoint_list = []                              #输出对应的element（调试用）

            for i in number_crack_ip:
                sdv_ele_list.append(sdv_element[i].data)
                # sdv_ele_lab_list.append(i.elementLabel)           #输出对应的element（调试用）
                # sdv_ele_ipoint_list.append(i.integrationPoint)    #输出对应的element（调试用）


            #sdv_element_lab1 = sdv_ele_lab_list [0::3]             #输出对应的element（调试用）
            sdv_ipoint_1  = sdv_ele_list [0::3]         #将sdv_element每3个元素中的第一个取出构成新数组，
            sdv_ipoint_2  = sdv_ele_list [1::3]
            sdv_ipoint_3  = sdv_ele_list [2::3]

            for i in range( 0,len(sdv_ipoint_1) ):
                status = sdv_ipoint_1[i] + sdv_ipoint_2[i] + sdv_ipoint_3[i]
                if status < 3:                   #当status小于设置的阀值时（checkValues=3），裂纹体积累加
                    crack_length  += evol_ele_list[i]*( 3-status )/3
            #print frameTime,crack_length                                   #提示信息面板中打印结果
            frame_time_list.append(frameTime)
            crack_length_list.append(crack_length)
            crack_data.append((frameTime,crack_length))
    totalTime += frame.frameValue                                           #step起始时间累加


    #计算不同frame裂纹面积的变化量
    delta_length_data=[]
    for i in range (len(crack_length_list)):
        try:
            delta_length  = crack_length_list[i+1] - crack_length_list[i]
            delta_length_data.append((frame_time_list[i+1],delta_length ))
        except:
            pass
    #----------------------------------------------------
    #endregion
    

    #region按CSV格式输出数据
    # f = file(output_path,'w')
    # f.write('step-time,crack_length\n')
    # for i in range(len(frame_time_list)) :
    #     data_line = '%f,%f\n'%(frame_time_list[i],crack_length_list[i])
    #     f.write(data_line)
    # f.close()

    # print 'The path of the result file is',output_path
    #endregion

    #region 处理数据输出到XYDATA
    #----------------------------------------------------
    import visualization

    unit_stime  = visualization.QuantityType(type=STEP_TIME)                    # 设置变量单位
    unit_arer   = visualization.QuantityType(type=AREA)
    unit_disp   = visualization.QuantityType(type=DISPLACEMENT)
    unit_volume = visualization.QuantityType(type=VOLUME)
    #创建xydata

    #---------------CLEN随时间的变化-----------
    cxy_name ='CLEN-time-'+ xy_name
    xyData = session.XYData(name= cxy_name, data=crack_data, sourceDescription='get from fieldOutputs', 
        contentDescription='total crack area ', positionDescription='for all cohesive element',
        axis1QuantityType = unit_stime, axis2QuantityType = unit_arer)

    #-----------CLEN的变化量随时间的变化--------
    dxy_name ='CLEN-D-time-'+ xy_name
    xyData = session.XYData(name= dxy_name, data=delta_length_data, sourceDescription='get from fieldOutputs', 
        contentDescription='total crack area ', positionDescription='for all cohesive element',
        axis1QuantityType = unit_stime, axis2QuantityType = unit_arer)


    ##写入历史变量
    # h_point  = HistoryPoint(coh_set)                                                                    #设置集合
    # h_region = step1.HistoryRegion(name='coh_set',description='all cohesive element',point=h_point)     #添加描述
    # h_cracklen = h_region.HistoryOutput(name='CLEN',description='Total area in coh_set',type=SCALAR)    #创建历史输出
    # h_cracklen.addData(frameValue = frame_time_list, value = crack_length_list)                         #向历史输出中添加数据
    # print 'HistoryOutput CLEN has been created'

    
    #处理xydata
    #----------CFN2-CUTTRT1---------
    try:
        xy_result = session.XYDataFromHistory(name='C1', 
            odb=odb, 
            outputVariableName='Total force due to contact pressure: CFN2     ASSEMBLY_CUTTER-1_SURF-M/ASSEMBLY_SET-CONTACT-C1_CNS_', 
            steps=('Step-load', ),
            useStepTime=True            )
    except :
        xy_result = session.XYDataFromHistory(name='C1', 
            odb=odb, 
            outputVariableName='Total force due to contact pressure: CFN2     ASSEMBLY_CUTTER-1_SURF-M/ASSEMBLY_PART-1-1_SET-G_CONTACT-C1_CNS_', 
            steps=('Step-load', ),
            useStepTime=True            )
    #----------CFN2-CUTTRT2---------
    try:
        xy_result = session.XYDataFromHistory(name='C2', 
            odb=odb, 
            outputVariableName='Total force due to contact pressure: CFN2     ASSEMBLY_CUTTER-2_SURF-M/ASSEMBLY_SET-CONTACT-C2_CNS_', 
            steps=('Step-load', ),
            useStepTime=True            )
    except :
        xy_result = session.XYDataFromHistory(name='C2', 
            odb=odb, 
            outputVariableName='Total force due to contact pressure: CFN2     ASSEMBLY_CUTTER-2_SURF-M/ASSEMBLY_PART-1-1_SET-G_CONTACT-C2_CNS_', 
            steps=('Step-load', ),
            useStepTime=True            )
    #----------U2-CUTTRT---------
    try:
        xy_result = session.XYDataFromHistory(name='U', 
        odb=odb, 
        outputVariableName='Spatial displacement: U2 PI: CUTTER-2 Node 9', 
        steps=('Step-load', ), 
        useStepTime=True            )
    except :
        xy_result = session.XYDataFromHistory(name='U', 
        odb=odb, 
        outputVariableName='Spatial displacement: U2 PI: CUTTER-2 Node 2', 
        steps=('Step-load', ), 
        useStepTime=True            ) 


    xy1 = session.xyDataObjects['C1']
    xy2 = session.xyDataObjects['C2']
    xy3 = session.xyDataObjects['U']
    xy4 = session.xyDataObjects[cxy_name]
    xy5 = session.xyDataObjects[dxy_name]

    #------Cutter_1-FN-----------
    xyc1 = combine(-xy3, xy1)
    xyc1.setValues(sourceDescription='Penetration depth vs Indentation force')
    tmpName = xyc1.name
    tmpName1 = 'Cutter_1-FN-'+ xy_name
    try:
        del session.xyDataObjects[tmpName1]
    except:
        pass
    session.xyDataObjects.changeKey(tmpName, tmpName1)

    #------Cutter_2-FN-----------
    xyc2 = combine(-xy3, xy2)
    xyc2.setValues(sourceDescription='Penetration depth vs Indentation force')
    tmpName = xyc2.name
    tmpName1 = 'Cutter_2-FN-'+ xy_name
    try:
        del session.xyDataObjects[tmpName1]
    except:
        pass 
    session.xyDataObjects.changeKey(tmpName, tmpName1) 

    #------------CLEN------------
    xyc3 = combine(-xy3, xy4)
    xyc3.setValues(sourceDescription='Penetration depth vs crack arae')
    tmpName = xyc3.name
    tmpName1 = 'CLEN-'+ xy_name
    try:
        del session.xyDataObjects[tmpName1]
    except:
        pass 
    session.xyDataObjects.changeKey(tmpName, tmpName1)

    # #------------CLEN-D-----------
    # xyc4 = combine(-xy3, xy5)
    # xyc4.setValues(sourceDescription='Penetration depth vs crack arae')
    # tmpName = xyc4.name
    # tmpName1 = 'CLEN-D-'+ xy_name
    # try:
    #     del session.xyDataObjects[tmpName1]
    # except:
    #     pass 
    # session.xyDataObjects.changeKey(tmpName, tmpName1)

    # session.xyDataObjects[dxy_name].setValues(axis1QuantityType=unit_disp, 
    #     axis2QuantityType=unit_volume)
    
    #---------Total-force---------
    xyc5 = xy1+xy2
    xyc5.setValues(sourceDescription='step time vs Indentation force')
    tmpName = xyc5.name
    tmpName1 = 'Total-force-time-'+ xy_name
    try:
        del session.xyDataObjects[tmpName1]
    except:
        pass 
    session.xyDataObjects.changeKey(tmpName, tmpName1)

    xyc6 = combine(-xy3, xy1+xy2)
    xyc6.setValues(sourceDescription='Penetration depth vs Indentation force')
    tmpName = xyc6.name
    tmpName1 = 'Total-force-'+ xy_name
    try:
        del session.xyDataObjects[tmpName1]
    except:
        pass 
    session.xyDataObjects.changeKey(tmpName, tmpName1) 


    del session.xyDataObjects['C1']
    del session.xyDataObjects['C2']
    del session.xyDataObjects['U']

    print '\nGet data from',odb_name
    print "\n==================Finshed=================="
    #----------------------------------------------------
    #endregion

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
    返回object\n
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
    fields = (('Multiselect mode',''),('Beginninig of slice',''),('End of slice',''))

    list_odject,slice_a,slice_b = getInputs ( fields=fields,label=title )


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

    return s_object
#================================================输入===============================================

output_path = 'G:\Desktop\output.csv' 
instName    = 'PART-1-1'                                        #岩石部件名称
setName     = 'COH3D6'                                          #包含所有cohesive单元的set名称
input_parameter =(output_path,instName,setName)


#=====main=====
def main ():
    t_start = time.time()
    
    obj  = session.odbs.values()
    name = 'Odbs'
    odb_list = Select_object(name,obj)

    for  i in odb_list:
        get_crack(i,input_parameter)

    t_end = time.time()
    print "Total time: %.2fs" %(t_end-t_start)
#==============

if __name__ =='__main__':
    print   '\n//////////作为主程序调用//////////'
    main()  # 作为主程序调用时执行
else:
    print '\n\n---------------------%s被其他程序调用---------------------'%__file__.split("\\")[-1]

