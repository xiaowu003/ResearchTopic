#!/user/bin/python
# -* - coding:UTF-8 -*- 
import os
from abaqusGui import getAFXApp
toolset=getAFXApp().getAFXMainWindow().getPluginToolset()

absPath = os.path.abspath(__file__)
absDir  = os.path.dirname(absPath)                      #获取当前文件的路径
#helpUrl = os.path.join(absDir, 'triad-help.txt')        #指定帮助文档

# Register a kernel plug-in in the Plug-ins menu.
#
toolset.registerKernelMenuButton(
    moduleName='crack_length', 
    functionName='main()',
    buttonText=('COH Crack|Statistical crack area'),
    applicableModules=['Visualization'],
    version='2.0', 
    author='Liu Xin',
    description='统计cohesive裂纹面积\n\n' \
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
    description='显示完全破坏的cohsive单元 \n\n' \
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
    description='批量提交job \n\n' \
                'a kernel command from a plug-in.\n' \
                "This plug-in's files may be copied from " + absDir,
    #helpUrl=helpUrl
)