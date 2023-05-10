---
title: 通过python读取EIGENVAL获取wannier90构造tb模型的冻窗口和解纠缠窗口
tags: Topology Python vasp
layout: article
license: true
toc: true
key: a20211222
pageview: true
# cover: /assets/images/GroupTheory/cube_symmetry.jpg
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
在通过Wannier90构建tb模型的时候，总是需要比较好的确定冻窗口和解纠缠窗口，而且这两个窗口的选取会很大的影响构造的tb的好坏，
这里就想通过代码，先粗略的估计一下这两个能量窗口的范围，再通过进一步的细致调节来构造处比较好的tb模型。
{:.info}
<!--more-->
# 获取价带和导带
在构造Wannier的时候，首先需要明确的就是想要构造那一部分能带的Wannier，通常就是费米面附近的能带，这对研究很多性质都是非常重要的。
首先是通过fatband分析，看看在费米面附近主要是哪些原子的那些轨道产生了贡献，以Bi为例，原胞中含有两个Bi原子，通过fatbang分析发现
费米面附近主要是由$s,p_x,p_y,p_z$轨道贡献，而且体系是考虑SOC的，所有在构建Wannier的时候需要考虑的能带的数目就是

$$
2*(1+1+1+1)*2 = 16
$$

这里的第一个2表示原胞中的两个Bi原子，四个1则是代表4条轨道，最后的2是因为考虑的SOC，所以还是需要考虑每条轨道的自旋自由度。

在确定了要构建的Wannier能带的数目之后，下面的一步就是要确定费米面的位置，这个值通过是可以在VASP计算的OUTCAR中得到的
```shell
cat OUTCAR |grep E-fermi|awk 'BEGIN {FS=" "} {print $3}' # 从OUTCAR中提取自洽计算的费米能
```
在得到了费米能之后，就可以进一步通过`EIGENVAL`中的结果来确定每个k点上的能量，将所有的k点连接起来就是一条能带，可以获取每条能带的带顶和带底，通过与费米能量的比较
就可以知道哪些带是占据能带(价带)，哪些带是空带(导带)。
```python
import numpy as np
import matplotlib.pyplot as plt
import os
def valread():
    """通过读取EIGENVAL文件来获取每条能带的最大值和最小值"""
    os.chdir(os.getcwd()) # 更改工作路径到当前程序所在路径
#     fermi = float(os.popen("sh /home/liyuxuan/script/fermi.sh").read()) # 运行提取费米能的脚本(在Linux下使用)
    fermi = 5.2759 # 自洽计算得到的费米能
    # 因为此时计算的时候是以费米面为基础的，所以要想得将化学势自然放在零，需要在能带计算中将费米能减去
    file = open('EIGENVAL') 
    temp = file.readline() # 先读取几个不太重要的行
    temp = file.readline() 
    temp = file.readline() 
    temp = file.readline() 
    system = file.readline() # INCAR中设置的SYSTEM
    info = file.readline() 
    info = [int(x) for x in info.strip().split()] # 第六行数据，第一个表示系统价电子数目，第二个表示k点数目，第三个表示能带数

    val = np.zeros((info[1],info[2]),dtype = float) # 构建数组来存储每个k点上每条能带的本征值
    klist = [] # 存储每个k点的值

    for i0 in range(info[1]): # 读取所有的k点
        f7 = file.readline()  # 读取一个空白行
        ktemp = file.readline()
        ktemp = [float(x) for x in ktemp.strip().split()] # 将字符串转换为浮点数
        klist.append(ktemp)
        for i1 in range(info[2]):
            temp1 = file.readline()
            temp1 = [float(x) for x in temp1.strip().split()]
            val[i0,i1] = temp1[1] # 只存储能带本征值
    # 对所有k点下，每条能带的范围进行分辨
    ic = 0
    iv = 0
    bandRange = np.zeros((info[2],2),dtype = float) # 存储每条能带的带顶和带底
    for i0 in range(info[2]): # 遍历每条能带
        temp_band = val[:,i0] - fermi
        band_min = np.min(temp_band)
        band_max = np.max(temp_band)
        bandRange[i0,0] = band_min #带底(相对而言，不是绝对的)
        bandRange[i0,1] = band_max #带顶
        if band_min > fermi:
            ic += 1
            print(" The energy of  %d(%d) th conduct band range is: %0.3f ---> %0.3f"%(ic,i0 + 1,band_min,band_max))
        else:
            iv += 1
            print(" The energy of  %d(%d) th valance band range is: %0.3f ---> %0.3f"%(iv,i0 + 1,band_min,band_max))

    return bandRange,val  #返回每条能带的范围和对应的的所有能带本征值的信息
```
通过上面的代码就可以自动的找到哪些带是占据态，哪些带是空态，这里在做的时候为了将费米面设置在零的位置，所以对每个k点的值都剪去了费米能。而费米能量的获取有两种方式，可以通过手动设置，也可以通过脚本自动实现
```shell
fermi = float(os.popen("sh path/fermi.sh").read()) # 运行提取费米能的脚本(在Linux下使用)
```
这里`fermi.sh`的内容如下
```shell
cat OUTCAR |grep E-fermi|awk 'BEGIN {FS=" "} {print $3}' # 从OUTCAR中提取自洽计算的费米能
```
可以通过python执行Linux的命令自动来获取，而得到了价带和导带之后，就可以来设置投影的时候需要使用的能量范围了。

# 确定能量窗口
通过前面的方式首先已经确定了所有的导带和价带，接下来就是根据前面计算得到的需要使用的Wannier函数的数目，来选择合适数量的导带和价带，这里需要16个Wannier函数，就直接讲它均分到导带和价带上面，所以可以选择8条导带和8条价带，而为了解纠缠可以选择分别预留
两条导带和两条价带。所以这里的思路就是以费米能量为基准，首先确定价带顶和导带底，分别向上、向下选择相同数目的导带和价带，以这个选择为基准，价带的带底作为冻窗口能带的最小值，导带的带顶作为冻窗口能量的最大值。至于解纠缠窗口的选择和冻窗口的选择方法是相同
的，唯一不同的也就是导带和价带选择数目的不同。
```python
def windows():
    """
    以费米能为基准，通过给定的能带数目来确定这个能带数所对应的整体能量范围
    """
    os.chdir(os.getcwd()) # 更改工作路径到当前程序所在路径
#     fermi = float(os.popen("sh /home/liyuxuan/script/fermi.sh").read()) # 运行提取费米能的脚本(在Linux下使用)
    fermi = 5.2759 # 自洽计算得到的费米能
    # 因为此时计算的时候是以费米面为基础的，所以要想得将化学势自然放在零，需要在能带计算中将费米能减去
    file = open('EIGENVAL') 
    temp = file.readline() # 先读取几个不太重要的行
    temp = file.readline() 
    temp = file.readline() 
    temp = file.readline() 
    system = file.readline() # INCAR中设置的SYSTEM
    info = file.readline() 
    info = [int(x) for x in info.strip().split()] # 第六行数据，第一个表示系统价电子数目，第二个表示k点数目，第三个表示能带数

    val = np.zeros((info[1],info[2]),dtype = float) # 构建数组来存储每个k点上每条能带的本征值
    klist = [] # 存储每个k点的值

    for i0 in range(info[1]): # 读取所有的k点
        f7 = file.readline()  # 读取一个空白行
        ktemp = file.readline()
        ktemp = [float(x) for x in ktemp.strip().split()] # 将字符串转换为浮点数
        klist.append(ktemp)
        for i1 in range(info[2]):
            temp1 = file.readline()
            temp1 = [float(x) for x in temp1.strip().split()]
            val[i0,i1] = temp1[1] # 只存储能带本征值
    # 对所有k点下，每条能带的范围进行分辨
    ic = 0
    iv = 0
    bandRange = np.zeros((info[2],2),dtype = float) # 存储每条能带的带顶和带底
    for i0 in range(info[2]): # 遍历每条能带，并输出所有能带的能量范围
        temp_band = val[:,i0] - fermi # VASP计算出来是存在化学势的，这里减去化学势就是为了将化学势放在零点
        band_min = np.min(temp_band)
        band_max = np.max(temp_band)
        bandRange[i0,0] = band_min #带底(相对而言，不是绝对的)
        bandRange[i0,1] = band_max #带顶
        if band_min > fermi:
            ic += 1
            print(" The energy of  %d(%d) th conduct band range is: %0.3f ---> %0.3f"%(ic,i0 + 1,band_min,band_max))
        else:
            iv += 1
            print(" The energy of  %d(%d) th valance band range is: %0.3f ---> %0.3f"%(iv,i0 + 1,band_min,band_max))
    
    # 确定能隙位于那两条能带之间
    fermi_band = 0 # 用这个参数来记录导带和价带分界时候的能带的指标
    for i0 in range(info[2]): # 遍历每条能带
        temp_band = val[:,i0] - fermi
        band_min = np.min(temp_band)
        band_max = np.max(temp_band)
        bandRange[i0,0] = band_min #带底(相对而言，不是绝对的)
        bandRange[i0,1] = band_max #带顶
        if band_min > fermi:
            fermi_band = i0 + 1 # 找到费米能所处的能带位置
            print("Fermi energy lying between %d th ------> %d th band" % (fermi_band - 1,fermi_band))
            break # 找到之后退出for循环

    # 根据给定的能带数，确定能带范围(冻窗口)
    frozeband = 7 # 选择加大和导带分别占据几条能带
    #engmin = bandRange[fermi_band - frozeband*2,0] + fermi # 这里又把化学势补了进去，对应的就是真实的VASP计算能带
    #engmax = bandRange[fermi_band + frozeband*2,1] + fermi
    engmin = bandRange[fermi_band - frozeband*2,0]
    engmax = bandRange[fermi_band + frozeband*2,1]
    print("Number of occ or unocc is %d ,energy range of froze windowns is : " % frozeband)
    print("%0.5f----->%0.5f" % (engmin,engmax))

    # 根据给定的能带数，确定解纠缠的能量窗口
    entband = frozeband + 1 # 根据价带和导带占据的数目，额外都多给一条作为解纠缠使用
    #engmin = bandRange[fermi_band - entband*2,0] + fermi # 这里又把化学势补了进去，对应的就是真实的VASP计算能带
    #engmax = bandRange[fermi_band + entband*2,1] + fermi
    engmin = bandRange[fermi_band - entband*2,0] # 这里又把化学势补了进去，对应的就是真实的VASP计算能带
    engmax = bandRange[fermi_band + entband*2,1]
    print("Number of occ or unocc is %d ,energy range of ent windows is : " % entband)
    print("%0.5f----->%0.5f" % (engmin,engmax))
```
在这里为了和前面分别哪些是导带哪些是价带的结果相同，就仍然是在数据中减去了费米能量减去了，这是因为我在做fatband分析的时候，就是以费米能量为基准进行的，所以这里也为了一致，就将费米能量减去了，通过这个简单的程序，就可以初步
估计在投影轨道的时候所需要设置的冻窗口的能量范围和解纠缠的能量范围，而为了让Wannier能带雨VASP能带你和的更好，还是需要进行一些微调。
# 完整代码
```python
import numpy as np
import matplotlib.pyplot as plt
import os
def valread():
    """通过读取EIGENVAL文件来获取每条能带的最大值和最小值"""
    os.chdir(os.getcwd()) # 更改工作路径到当前程序所在路径
#     fermi = float(os.popen("sh /home/liyuxuan/script/fermi.sh").read()) # 运行提取费米能的脚本(在Linux下使用)
    fermi = 5.2759 # 自洽计算得到的费米能
    # 因为此时计算的时候是以费米面为基础的，所以要想得将化学势自然放在零，需要在能带计算中将费米能减去
    file = open('EIGENVAL') 
    temp = file.readline() # 先读取几个不太重要的行
    temp = file.readline() 
    temp = file.readline() 
    temp = file.readline() 
    system = file.readline() # INCAR中设置的SYSTEM
    info = file.readline() 
    info = [int(x) for x in info.strip().split()] # 第六行数据，第一个表示系统价电子数目，第二个表示k点数目，第三个表示能带数

    val = np.zeros((info[1],info[2]),dtype = float) # 构建数组来存储每个k点上每条能带的本征值
    klist = [] # 存储每个k点的值

    for i0 in range(info[1]): # 读取所有的k点
        f7 = file.readline()  # 读取一个空白行
        ktemp = file.readline()
        ktemp = [float(x) for x in ktemp.strip().split()] # 将字符串转换为浮点数
        klist.append(ktemp)
        for i1 in range(info[2]):
            temp1 = file.readline()
            temp1 = [float(x) for x in temp1.strip().split()]
            val[i0,i1] = temp1[1] # 只存储能带本征值
    # 对所有k点下，每条能带的范围进行分辨
    ic = 0
    iv = 0
    bandRange = np.zeros((info[2],2),dtype = float) # 存储每条能带的带顶和带底
    for i0 in range(info[2]): # 遍历每条能带
        temp_band = val[:,i0] - fermi
        band_min = np.min(temp_band)
        band_max = np.max(temp_band)
        bandRange[i0,0] = band_min #带底(相对而言，不是绝对的)
        bandRange[i0,1] = band_max #带顶
        if band_min > fermi:
            ic += 1
            print(" The energy of  %d(%d) th conduct band range is: %0.3f ---> %0.3f"%(ic,i0 + 1,band_min,band_max))
        else:
            iv += 1
            print(" The energy of  %d(%d) th valance band range is: %0.3f ---> %0.3f"%(iv,i0 + 1,band_min,band_max))

    return bandRange,val  #返回每条能带的范围和对应的的所有能带本征值的信息
#----------------------------------------------------------------------------------------------------------------
def windows():
    """
    以费米能为基准，通过给定的能带数目来确定这个能带数所对应的整体能量范围
    """
    os.chdir(os.getcwd()) # 更改工作路径到当前程序所在路径
#     fermi = float(os.popen("sh /home/liyuxuan/script/fermi.sh").read()) # 运行提取费米能的脚本(在Linux下使用)
    fermi = 5.2759 # 自洽计算得到的费米能
    # 因为此时计算的时候是以费米面为基础的，所以要想得将化学势自然放在零，需要在能带计算中将费米能减去
    file = open('EIGENVAL') 
    temp = file.readline() # 先读取几个不太重要的行
    temp = file.readline() 
    temp = file.readline() 
    temp = file.readline() 
    system = file.readline() # INCAR中设置的SYSTEM
    info = file.readline() 
    info = [int(x) for x in info.strip().split()] # 第六行数据，第一个表示系统价电子数目，第二个表示k点数目，第三个表示能带数

    val = np.zeros((info[1],info[2]),dtype = float) # 构建数组来存储每个k点上每条能带的本征值
    klist = [] # 存储每个k点的值

    for i0 in range(info[1]): # 读取所有的k点
        f7 = file.readline()  # 读取一个空白行
        ktemp = file.readline()
        ktemp = [float(x) for x in ktemp.strip().split()] # 将字符串转换为浮点数
        klist.append(ktemp)
        for i1 in range(info[2]):
            temp1 = file.readline()
            temp1 = [float(x) for x in temp1.strip().split()]
            val[i0,i1] = temp1[1] # 只存储能带本征值
    # 对所有k点下，每条能带的范围进行分辨
    ic = 0
    iv = 0
    bandRange = np.zeros((info[2],2),dtype = float) # 存储每条能带的带顶和带底
    for i0 in range(info[2]): # 遍历每条能带，并输出所有能带的能量范围
        temp_band = val[:,i0] - fermi # VASP计算出来是存在化学势的，这里减去化学势就是为了将化学势放在零点
        band_min = np.min(temp_band)
        band_max = np.max(temp_band)
        bandRange[i0,0] = band_min #带底(相对而言，不是绝对的)
        bandRange[i0,1] = band_max #带顶
        if band_min > fermi:
            ic += 1
            print(" The energy of  %d(%d) th conduct band range is: %0.3f ---> %0.3f"%(ic,i0 + 1,band_min,band_max))
        else:
            iv += 1
            print(" The energy of  %d(%d) th valance band range is: %0.3f ---> %0.3f"%(iv,i0 + 1,band_min,band_max))
    
    # 确定能隙位于那两条能带之间
    fermi_band = 0 # 用这个参数来记录导带和价带分界时候的能带的指标
    for i0 in range(info[2]): # 遍历每条能带
        temp_band = val[:,i0] - fermi
        band_min = np.min(temp_band)
        band_max = np.max(temp_band)
        bandRange[i0,0] = band_min #带底(相对而言，不是绝对的)
        bandRange[i0,1] = band_max #带顶
        if band_min > fermi:
            fermi_band = i0 + 1 # 找到费米能所处的能带位置
            print("Fermi energy lying between %d th ------> %d th band" % (fermi_band - 1,fermi_band))
            break # 找到之后退出for循环

    # 根据给定的能带数，确定能带范围(冻窗口)
    frozeband = 7 # 选择加大和导带分别占据几条能带
    #engmin = bandRange[fermi_band - frozeband*2,0] + fermi # 这里又把化学势补了进去，对应的就是真实的VASP计算能带
    #engmax = bandRange[fermi_band + frozeband*2,1] + fermi
    engmin = bandRange[fermi_band - frozeband*2,0]
    engmax = bandRange[fermi_band + frozeband*2,1]
    print("Number of occ or unocc is %d ,energy range of froze windowns is : " % frozeband)
    print("%0.5f----->%0.5f" % (engmin,engmax))

    # 根据给定的能带数，确定解纠缠的能量窗口
    entband = frozeband + 1 # 根据价带和导带占据的数目，额外都多给一条作为解纠缠使用
    #engmin = bandRange[fermi_band - entband*2,0] + fermi # 这里又把化学势补了进去，对应的就是真实的VASP计算能带
    #engmax = bandRange[fermi_band + entband*2,1] + fermi
    engmin = bandRange[fermi_band - entband*2,0] # 这里又把化学势补了进去，对应的就是真实的VASP计算能带
    engmax = bandRange[fermi_band + entband*2,1]
    print("Number of occ or unocc is %d ,energy range of ent windows is : " % entband)
    print("%0.5f----->%0.5f" % (engmin,engmax))
#----------------------------------------------------------------------------------------------------------------
def main():
    # bandRange,val = valread() # 首先得到每条能带的取值范围和所有能带的信息
    windows()

#------------------------------------------------------------------------------------------------------------------
if __name__=="__main__":
    main()
```
# 参考
- 1.[Wannier90输入文件中num_wann, num_bands, 和energy window等参数设置规则
](https://wap.sciencenet.cn/blog-2909108-1154273.html?mobile=1)


# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)