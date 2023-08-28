#!/user/bin/python
# -*- coding: utf-8 -*-
#liuxin66@live.com
# ͳ��coh�������
# 2020��7��16�� 
# ���з�ʽ��
#   1���˵���File->Run Script...
#================================================================================
from abaqus import *
from abaqusConstants import *
from odbAccess import *
import time
#==============================================��������==============================================
def get_crack(odb_num,input_parameter):
    output_path = input_parameter[0] 
    instName    = input_parameter[1]                                       
    setName     = input_parameter[2] 

    odb = session.odbs.values()[odb_num]                            #ѡ���Ѵ򿪵�odb�еĵ�һ��
    myAssembly = odb.rootAssembly
    inst    = myAssembly.instances[instName]
    coh_set = inst.elementSets[setName]                             #ָ��set
    step1   = odb.steps.values()[0]                                 #ָ��������
    print "\n==================started=================="
    odb_name = odb.name.split('/')[-1].split('.')[0]
    l_name  = odb_name.split('-')
    xy_name = l_name[1]+'-'+l_name[-1]                               #xydata ���ݺ�׺
    print xy_name
    print 'Get data from',odb_name

    #region ��ȡ���һ�������������������ȷ�����Ƶ�ͳ������(ֻ�Ի���ѵ�coh��Ԫ����ͳ�����Ч��)
    #----------------------------------------------------
    final_frame_sdv = step1.frames[-1].fieldOutputs['SDV1'].getSubset(region = coh_set).values
    final_sdv_list  = []
    final_ele_lab_list = []
    for i in final_frame_sdv:                   #��sdvֵ��˳��洢���б��У����ڴ���
        final_sdv_list.append(i.data)
        #final_ele_lab_list.append(i.elementLabel)

    f_sdv_ip_1  = final_sdv_list [0::3]         #��sdv_elementÿ3��Ԫ���еĵ�һ��ȡ�����������飨�����ֵ�ֳ�3���б�
    f_sdv_ip_2  = final_sdv_list [1::3]
    f_sdv_ip_3  = final_sdv_list [2::3]
    #f_ele_lab   = final_ele_lab_list[0::3]
    max_crack   = len(f_sdv_ip_1)


    number_crack_ip   = []
    number_crack_ele  = []                      #ͳ������������ evol�б�� ���
    #lab_crack_ele     = []                      #ͳ������������ evol�б�� ��Ԫ���
    for i in range( 0,max_crack ):              #��������coh��Ԫ
        status = f_sdv_ip_1[i] + f_sdv_ip_2[i] + f_sdv_ip_3[i]   

        if status < 3:                          #���е�Ԫ�����ƻ�ʱ

            number_crack_ele.append(i)          #ͳ������������ evol�б�� ���
            #lab_crack_ele.append(f_ele_lab[i])  #ͳ������������ evol�б�� ��Ԫ���

            n=3*i             
            number_crack_ip.extend([n,n+1,n+2]) #ͳ������������ sdv1�б�� ���
            
    print 'MAX %s cohsive is cracked ' %max_crack
    #----------------------------------------------------
    #endregion

    #region ѭ����ȡ��ͬʱ�̵����
    #----------------------------------------------------
    totalTime = 0.0                                                 #��¼������ʱ��
    count=0                                                         #��¼����
    frame_time_list=[]
    crack_length_list=[]
    crack_data=[]
    evol_element = step1.frames[0].fieldOutputs['EVOL'].getSubset(region = coh_set).values
    evol_ele_list =[]
    # evol_ele_lab_list =[]                                         #�����Ӧ��element�������ã�
    for i in number_crack_ele:                                      #������ת����list
        evol_ele_list.append(evol_element[i].data)
        # evol_ele_lab_list.append(i.elementLabel)

    #frame= step1.frames[2]
    #ָ��evol��֡������֡�ڶ�һ��������ֻ��һ�����Ч�ʣ�
    #for step in odb.steps.values():                                #�������з�����
    for frame in step1.frames:                                      #��������֡
            print step1.name,'frame=',count,'...'
            count += 1
            crack_length = 0                                        #ÿ֡���Ƴ��ȸ���ֵ
            frameTime   = frame.frameValue + totalTime              #��ǰ֡��Ӧ����ʱ��
            
            sdv_element   = frame.fieldOutputs['SDV1'].getSubset(region = coh_set).values #sdv������������ǻ��ֵ��ϵ�ֵ��һ��coh3d6��3�����ֵ�

            #������ת����list
            sdv_ele_list        = []
            # sdv_ele_lab_list    = []                              #�����Ӧ��element�������ã�
            # sdv_ele_ipoint_list = []                              #�����Ӧ��element�������ã�

            for i in number_crack_ip:
                sdv_ele_list.append(sdv_element[i].data)
                # sdv_ele_lab_list.append(i.elementLabel)           #�����Ӧ��element�������ã�
                # sdv_ele_ipoint_list.append(i.integrationPoint)    #�����Ӧ��element�������ã�


            #sdv_element_lab1 = sdv_ele_lab_list [0::3]             #�����Ӧ��element�������ã�
            sdv_ipoint_1  = sdv_ele_list [0::3]         #��sdv_elementÿ3��Ԫ���еĵ�һ��ȡ�����������飬
            sdv_ipoint_2  = sdv_ele_list [1::3]
            sdv_ipoint_3  = sdv_ele_list [2::3]

            for i in range( 0,len(sdv_ipoint_1) ):
                status = sdv_ipoint_1[i] + sdv_ipoint_2[i] + sdv_ipoint_3[i]
                if status < 3:                   #��statusС�����õķ�ֵʱ��checkValues=3������������ۼ�
                    crack_length  += evol_ele_list[i]*( 3-status )/3
            #print frameTime,crack_length                                   #��ʾ��Ϣ����д�ӡ���
            frame_time_list.append(frameTime)
            crack_length_list.append(crack_length)
            crack_data.append((frameTime,crack_length))
    totalTime += frame.frameValue                                           #step��ʼʱ���ۼ�


    #���㲻ͬframe��������ı仯��
    delta_length_data=[]
    for i in range (len(crack_length_list)):
        try:
            delta_length  = crack_length_list[i+1] - crack_length_list[i]
            delta_length_data.append((frame_time_list[i+1],delta_length ))
        except:
            pass
    #----------------------------------------------------
    #endregion
    

    #region��CSV��ʽ�������
    # f = file(output_path,'w')
    # f.write('step-time,crack_length\n')
    # for i in range(len(frame_time_list)) :
    #     data_line = '%f,%f\n'%(frame_time_list[i],crack_length_list[i])
    #     f.write(data_line)
    # f.close()

    # print 'The path of the result file is',output_path
    #endregion

    #region �������������XYDATA
    #----------------------------------------------------
    import visualization

    unit_stime  = visualization.QuantityType(type=STEP_TIME)                    # ���ñ�����λ
    unit_arer   = visualization.QuantityType(type=AREA)
    unit_disp   = visualization.QuantityType(type=DISPLACEMENT)
    unit_volume = visualization.QuantityType(type=VOLUME)
    #����xydata

    #---------------CLEN��ʱ��ı仯-----------
    cxy_name ='CLEN-time-'+ xy_name
    xyData = session.XYData(name= cxy_name, data=crack_data, sourceDescription='get from fieldOutputs', 
        contentDescription='total crack area ', positionDescription='for all cohesive element',
        axis1QuantityType = unit_stime, axis2QuantityType = unit_arer)

    #-----------CLEN�ı仯����ʱ��ı仯--------
    dxy_name ='CLEN-D-time-'+ xy_name
    xyData = session.XYData(name= dxy_name, data=delta_length_data, sourceDescription='get from fieldOutputs', 
        contentDescription='total crack area ', positionDescription='for all cohesive element',
        axis1QuantityType = unit_stime, axis2QuantityType = unit_arer)


    ##д����ʷ����
    # h_point  = HistoryPoint(coh_set)                                                                    #���ü���
    # h_region = step1.HistoryRegion(name='coh_set',description='all cohesive element',point=h_point)     #�������
    # h_cracklen = h_region.HistoryOutput(name='CLEN',description='Total area in coh_set',type=SCALAR)    #������ʷ���
    # h_cracklen.addData(frameValue = frame_time_list, value = crack_length_list)                         #����ʷ������������
    # print 'HistoryOutput CLEN has been created'

    
    #����xydata
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

#��ȡҪ�����ģ�ͷ�Χ
def simple (name):
    '''
    ��odbname��ȥ��·�����õ��ļ���
    '''
    if '/' in name :                                       
        simple_name = name.split('/')[-1].split('.')[0]
    else:
        simple_name = name
    return  simple_name
def Select_object(name,objects):
    '''
    ѡ��Ҫ����Ķ��󣬺ʹ���˳��
    ����object\n
    ʹ��ʾ����
           name ='JOBs'
           obj = mdb.jobs.values()
           (joblist,viewlist) = Select_object(name,obj)
           ������Ϊ����ȫѡ����object
    '''
    title   = 'Select the %s from below\n\n'  %name
    count   = 0
    #���job�б�
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


    #���object�����б�
    if list_odject == '':                               #���Ϊ����slice�����б�
        if slice_a == '' and slice_b == '':
            s_object = range(count)
        else:
            slice_a = int (slice_a)
            slice_b = int (slice_b) + 1

            s_object   = range(slice_a,slice_b,1)
    elif list_odject =='all':
        s_object = range(count)
    else:
        s_object = map(int,list_odject.split(','))      #object ��0��ʼ���

    return s_object
#================================================����===============================================

output_path = 'G:\Desktop\output.csv' 
instName    = 'PART-1-1'                                        #��ʯ��������
setName     = 'COH3D6'                                          #��������cohesive��Ԫ��set����
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
    print   '\n//////////��Ϊ���������//////////'
    main()  # ��Ϊ���������ʱִ��
else:
    print '\n\n---------------------%s�������������---------------------'%__file__.split("\\")[-1]

