#!/user/bin/python
# -* - coding:UTF-8 -*- 
import os
from abaqusGui import getAFXApp
toolset=getAFXApp().getAFXMainWindow().getPluginToolset()

absPath = os.path.abspath(__file__)
absDir  = os.path.dirname(absPath)                      #��ȡ��ǰ�ļ���·��
#helpUrl = os.path.join(absDir, 'triad-help.txt')        #ָ�������ĵ�

# Register a kernel plug-in in the Plug-ins menu.
#
toolset.registerKernelMenuButton(
    moduleName='crack_length', 
    functionName='main()',
    buttonText=('COH Crack|Statistical crack area'),
    applicableModules=['Visualization'],
    version='2.0', 
    author='Liu Xin',
    description='ͳ��cohesive�������\n\n' \
                 'a kernel command from a plug-in.\n' \
                 "This plug-in's files may be copied from " + absDir,
    #helpUrl=helpUrl
)

toolset.registerKernelMenuButton(
    moduleName='display_crack', 
    functionName='main()',
    buttonText=('COH Crack|Display_crack'),
    applicableModules=['Visualization'],
    version='2.0', 
    author='Liu Xin',
    description='��ʾ��ȫ�ƻ���cohsive��Ԫ \n\n' \
                'a kernel command from a plug-in.\n' \
                "This plug-in's files may be copied from " + absDir,
    #helpUrl=helpUrl
)

toolset.registerKernelMenuButton(
    moduleName='job_submit', 
    functionName='main()',
    buttonText=('job_submit'),
    applicableModules=['Job'],
    version='1.0', 
    author='Liu Xin',
    description='�����ύjob \n\n' \
                'a kernel command from a plug-in.\n' \
                "This plug-in's files may be copied from " + absDir,
    #helpUrl=helpUrl
)