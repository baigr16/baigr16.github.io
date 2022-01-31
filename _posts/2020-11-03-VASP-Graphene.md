---
title: VASP计算Graphene能带
tags: vasp
layout: article
license: true
toc: true
key: a20201103
pageview: true
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
最近看了很多和第一性原理计算相关的资料,知识总是要和实践结合,这里就利用VASP来计算一下Graphene的能带,作为研究很透彻的一个体系,起码在$\Gamma$点的Dirac锥的结构势清楚了,所以就拿这个体系来入门.下面的内容都是我最为一个初学者的理解,其中应该会包含一些错误的理解,希望可以指出.
{:.info}
<!--more-->
# 1.晶格弛豫
在对一个体系进行计算时,首先是要确定你的这个体系时稳定的,所以要先进行晶格弛豫,我理解的就是**晃动一下整个晶格,重要保证它有一定的稳定性,不至于轻微的晃动就变成另外一种结构.**这里这计算的时2维的Graphene,虽然从结构上看确实好像有一定稳定性,但是还是要进行弛豫计算.
## INCAR
```python
ISTART=0   # 参数0代表开始一个全新的计算
ICHARG=2   # 2代表初始电荷密度由原子轨道电荷密度叠加构成
LWAVE=.FALSE.  # 轨道以及电荷密度不进行输出(节省空间,但有时你可能需要这些数据,改为true即可)
LCHARG=.FALSE. # 控制电荷摸底输出

PREC=Normal# standard precision   计算精度控制
ENMAX=450# cut off set according to POTCAR   最大截断能
EDIFF=1E-5# electron step accurate   电子步计算时的精度控制
ISMEAR=0    # 占据函数设置,0为高斯占据函数
SIGMA=0.05  # 通常ismear为0时,使用sigma=0.05

IBRION=2# use GGA algorithm   利用共轭梯度法计算离子弛豫
ISIF=3# relax both cell and atom   力及力张量计算,3代表全部输出
NSW=100# 100 ionic step   最大离子步数目
EDIFFG=-0.03# forces smaller 0.02 A/eV   离子步计算精度
```
## POSCAR
关于POSCAR的设置可以参考这篇博客[VASP基本输入文件准备](https://yxli8023.github.io/2020/09/27/VASP-4.html)
```python
graphene
1.0
2.4677236   0.0000000   0.0000000
-1.2338618   2.1371113   0.0000000
0.0000000   0.0000000  15.0000000
C
2
Direct
0.00000000   0.00000000   0.25000000
0.66666667   0.33333333   0.25000000
```
## KPOINTS
```python
AUTO
0
Gamma
6 6 1
0 0 0
```

晶格弛豫过程主要是为了确定晶格的一个稳定结构,计算完成之后需要的只是CONTCAR文件,它包含了最终计算得到的稳定结果

![png](/assets/images/vasp/contcar.png)

这里第三个基矢发生了一些小的变化,但是这里计算的是2维Graphene,还是将它修改维最初的15.00000

# 自洽计算
完成了晶格弛豫稳定性计算之后,接下来就进行自洽计算,这个时候需要将上一步计算得到的CONTCAR复制为POSCAT来进行自洽计算
## INCAR
```python
ISTART = 0   # 还是要重新开始一个计算
ICHARG = 2   # 2代表初始电荷密度由原子轨道电荷密度叠加构成

PREC = Accurate    # accurate precision
ENMAX = 450
EDIFF = 1E-6
ISMEAR = 0
SIGMA = 0.05
               # 因为已经通过弛豫过程确定了稳定结构,所以这时候不再进行离子弛豫
IBRION = -1    # No ionic movement 
NSW = 0        # No ionic movement
```
## KPOINTS
```python
AUTO
0
Gamma    # 以Gamma点为中心自动生成k格点
12 12 1
0 0 0
```
自洽计算完成之后会得到电荷密度文件CHGCAR,利用这个文件来进行之后的能带计算.

# 能带计算
利用前两步计算得到的CHGCAR以及CONTCAR文件,将CONTCAR复制为POSCAR,进行能带计算.
## INCAR
```python
ICHARG = 11   # Read from CHGCAR  从现有的文件中读取电荷密度

PREC = Accurate  # accurate precision
ENMAX = 450   # cut off set according POTCAT
EDIFF = 1E-6
ISMEAR = 0
SIGMA = 0.05

IBRION = -1  # no ionic movement
NSW = 0      # no ionic movement
 
LORNIT = 11   # projected band   这里会计算每个轨道对能带的贡献
```


## KPOINTS
```python
K-points for bandstructure G-K-M-G
30
line-mode
reciprocal    # 下面就是设置高对称点的路径
0.000 0.000 0.000 # G
-0.333 0.666 0.000 # K

-0.333 0.666 0.000 # K
0.500 0.000 0.000 # M
 
0.500 0.000 0.000 # M
0.000 0.000 0.000 # G
```
执行完毕之后即可以得到最终的结果vasprun.xml,可以通过p4vasp来查看能带结果.

![png](/assets/images/vasp/graphene.png)

# 计算平台
我这里是在[并行云](https://www.paratera.com/companyNews/companyDetail18.html)的超算上进行了,自己提供VASP的包之后,平台会有人进行编译配置,完成之后会提供一个提交作业的脚本给你,展示一下我的
```shell
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
source /public3/soft/other/module.sh
module load mpi/intel/17.0.5-cjj-public3
export PATH=/public3/home/username/package/VASP/vasp.6.1.0/bin:$PATH
srun vasp_std
```

提交作业的时候,将这个脚本(比如名称为job.sh)复制到你的文件夹里面,然后执行

> sbatch job.sh

关于脚本与服务器的其他内容,如果取申请自然会有相关的资料给你.


# 参考
- 1.[vasp tutorial: 12. graphene band structure (with projection on atomic orbitals)](https://www.youtube.com/watch?v=I9uB-px_YUY&ab_channel=QuantumNerd)
