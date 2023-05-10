---
title: WannierTools研究Haldane 模型
tags:  vasp
layout: article
license: true
toc: true
key: a20210103
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
pageview: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
通过第一性原理计算,可以对具体的材料进行研究,我自己也在对这个方向慢慢进行摸索,但是VASP等软件并不是万能的,有些其它的性质它并不能得到,所以就有了一些开源的软件,比如[Wannier90](http://www.wannier.org/),[WannierTools](https://www.wanniertools.org/)可以用来你和紧束缚能带,计算体系的拓扑性质.我在这里首先学习的是WannierTools,因为它可以利用一些有Wannier90计算得到的数据,来计算对应体系的拓扑量以及能谱,自己对这方面是比较熟悉的,就先学习如何利用WannierTools来计算[Haldane模型](https://topocondmat.org/w4_haldane/haldane_model.html)的一些拓扑性质.
{:.info}
<!--more-->
首先如何下载和安装[WannierTools](https://www.wanniertools.org/)可以自行参考官网,通常情况下如果课题组从事这方面的研究肯定服务器上就会安装好,除非你是白手起家,或者是自己想学这个东西(比如我),那么你就需要去按照教程安装一下了,这里就不介绍怎么安装了,后面有机会我再整理一份教程来演示如何编译安装[WannierTools](https://www.wanniertools.org/).

这里的对与Haldane模型的研究,也是基于[WannierTools](https://www.wanniertools.org/)自带的练习进行的,我只是简单的跟着重复一些这些例子,然后加上一些自己对一些参数设置的理解.软件下载之后解压之后,主体结构如下

![png](/assets/images/wannierTools/f1.png)

文件夹`src`中放置的是计算的主体源代码,如果对某一部分内容感兴趣,可以直接研究学习源代码,`examples`中放置的就是一些具体的计算实例了,初学的话可以跟着这里面的例子来学习,官网也有这些例子的详细教程.

![png](/assets/images/wannierTools/f2.png)

# Haldane model

![png](/assets/images/wannierTools/Haldane1.png)

关于Haldane模型的参数和拓扑性质如上图所示,我这里主要就是利用`examples`中的`Haldane_model`这个文件夹中的内容来对这个实例进行重复,文件夹内容如下图所示

![png](/assets/images/wannierTools/f3.png)

在利用WannierTools计算的时候,主要的文件有两个,一个是`name_hr.dat`这个文件时候Wannier90计算产生的紧束缚模型的数据,当然了也可以通过现有的紧束缚模型来产生这个数据,实例中的`Haldane_model`就是利用现有的紧束缚模型来产生的这个数据,这个方法我还在研究中,至于如何利用Wannier90来产生这个数据,这个任务是之后的事情了.第二个重要的文件就是`wt.in`,这个文件的主要作用就是来进行计算参数的设置,控制要进行什么样的计算.
{:.warning}
首先进入到上图的文件夹中,执行下面的命令来产生上面提及到的两个重要文件
```python
python haldane_hr_gen-trivial-insulator.py
cp wt.in-trivial-insulator wt.in
```
如何利用紧束缚模型产生`name_hr.dat`我并没学会,所以这里主要就解读一下`wt.in`中的参数都是在做什么
```python
&TB_FILE
Hrfile = "Haldane_hr.dat"   ! 设置紧束缚近似模型数据存储的文件名
/


!> bulk band structure calculation flag
&CONTROL                  ! 从这里开始设置计算的细节
BulkBand_calc         = T  ! calculate band structure in kpath mode(计算体态能带)
BulkBand_plane_calc   = T  ! calculate band structure in kplane mode
SlabBand_calc         = T  ! Calculate slab band structure in kpath mode(半开边界的能带计算)
SlabSS_calc           = F  ! Calculate surface states in kpath mode(计算表面态,也就是计算表面的谱函数)
Wanniercenter_calc    = T  ! calculate Wilson loop(利用/wilson loop的方法计算体系Wannier Center的演化,可用来研究拓扑不变量)
BerryCurvature_calc   = T  ! 计算体系的Berry曲率
/

&SYSTEM
NSLAB =60               ! 半开边界情况下,开边界方向上格点的数目
NumOccupied = 1         ! NumOccupied  占据态的数目
SOC = 0                 ! soc(0代表没有自旋轨道耦合,1代表有自旋轨道后河)
E_FERMI = 0        ! e-fermi  费米能的位置
/

&PARAMETERS
Eta_Arc = 0.01     ! infinite small value, like brodening 
E_arc = 0.0         ! energy for calculate Fermi Arc
OmegaNum = 1000  ! omega number      设置确定能量区间中的撒点数目
OmegaMin = -5.0     ! energy interval  设置计算的能量区间
OmegaMax =  5.0     ! energy interval
Nk1 = 60            ! number k points 在k空间计算时候的k点撒点数目控制
Nk2 = 60            ! number k points 
NP = 1              ! number of principle layers
/

LATTICE  ! 这里设置元胞的信息,与VASP中的设置是相同的
Angstrom
2.1377110  -1.2342080   0.0000000
0.0000000   2.4684160   0.0000000
0.0000000   0.0000000   10.000000

ATOM_POSITIONS
2                               ! number of atoms for projectors
Direct                          ! Direct or Cartisen coordinate
C 0.333333 0.666667 0.500000    ! 元胞中原子的位置
C 0.666667 0.333333 0.500000 

PROJECTORS   ! 设置元胞中每个原子贡献的轨道
1 1          ! number of projectors
C pz
C pz


SURFACE            ! See doc for details(表面态计算设置,也就是选择哪个方向是开边界,哪个方向是周期的)
 0  0  1           ! 通常体系默认为是3D的,所以一般情况下是设置两个方向为周期,一个方向为开边界来计算能带
 1  0  0
 0  1  0

KPATH_BULK            ! k point path(这里用来设置计算体态能带时候的路径)
3              ! number of k line only for bulk band
  M   0.50000  0.00000  0.00000   K' -.33333   -.33333  0.00000
  K'  -.33333  -.33333  0.00000   G  0.00000   0.00000  0.00000
  G   0.00000  0.00000  0.00000   K  0.33333   0.33333  0.00000

KPATH_SLAB
1        ! numker of k line for 2D case
0 0.0 0.0 1 0. 1.0  ! k path for 2D case

KPLANE_SLAB
-0.5 -0.5      ! Original point for 2D k plane
 1.0  0.0      ! The first vector to define 2D k plane 
 0.0  1.0      ! The second vector to define 2D k plane  for arc plots

KPLANE_BULK
 0.00  0.00  0.00   ! Original point for 3D k plane 
 1.00  0.00  0.00   ! The first vector to define 3d k space plane
 0.00  1.00  0.00   ! The second vector to define 3d k space plane
```
将上面所说的两个文件准备好之后,就可以开始进行计算了
```python
mpirun -np 2 wt.x &  # -np后面的数值是你想要开多少个核进行并行计算
```
计算完成之后,根据设置的参数不同,则会有一些不同的文件结果输出,所有的结果都可以利用`gnuplot`来绘制,所以最后的结果中可以看到一些后缀为`.gnu`的文件,执行这个文件即可
```python
gnuplot bulkek.gnu
gnuplot bulkek_plane.gnu
gnuplot Berrycurvature.gnu
gnuplot wcc.gnu
gnuplot slabek.gnu
```
最后的结果如下所示
![png](/assets/images/wannierTools/ha2.png)

![png](/assets/images/wannierTools/ha3.png)

![png](/assets/images/wannierTools/ha4.png)

![png](/assets/images/wannierTools/ha5.png)

# Other Case
上面是对平庸相的Haldane模型计算的结果,在对参数进行调整之后就可以计算其它相对应的结果
- Chern insulator
```python
python haldane_hr_gen-chern-insulator.py
cp wt.in-chern-insulator wt.in
mpirun -np 2 wt.x 
```

- Gapless semimeta
```python
python haldane_hr_gen-gapless.py
cp wt.in-gapless wt.in
mpirun -np 2 wt.x &
```

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)