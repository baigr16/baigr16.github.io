---
title: wannierTools研究Topological Insulator
tags:  vasp
layout: article
license: true
toc: true
key: a20210103a
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
这里我想继续通过[WannierTools](https://www.wannierTools.org/)来计算Topological insulator的一些性质,来让自己对输入文件中的一些参数有一个更加深刻的认识,并通过这个实例来学习新的一些控制参数的作用.
{:.info}
<!--more-->
# Topological insulator
这个博客中学习的实例同样来自于源代码中`examples`文件夹中

![png](/assets/images/wannierTools/ti1.png)

![png](/assets/images/wannierTools/ti2.png)

首先解压文件得到这个体系的紧束缚模型的数据
```python
tar xzvf wannier90_hr.dat.tar.gz
```
`wt.in`已经准备好了,接下来就主要分析学习一下里面的参数设置

```fortran
&TB_FILE
Hrfile = 'wannier90_hr.dat'  ! 存储紧束缚能带数据的文件
Package = 'VASP'     
/


!> bulk band structure calculation flag
&CONTROL
BulkBand_calc         = T  ! 选择是否计算体态能带
SlabBand_calc         = F  ! Flag for 2D slab energy band calculation
SlabSS_calc           = T  ! Flag for surface state ARPES spectrum calculation
wanniercenter_calc    = T  ! 计算体系能带的Wilson loop,可用来判断体系拓扑性质
/

&SYSTEM
NSLAB = 10   ! Number of slabs for 2d Slab system
NumOccupied = 28        ! Number of occupied bands for bulk unit cell
SOC = 1                 ! A parameter to control soc;  Soc=0 means no spin-orbit coupling; Soc>0 means spin-orbit coupling
E_FERMI = -2.9600       ! Fermi energy, search E-fermi in OUTCAR for VASP, set to zero for Wien2k
/

&PARAMETERS
OmegaNum = 301     ! The number of energy slices between OmegaMin and OmegaMax      
OmegaMin = -1.0     ! energy interval
OmegaMax =  0.5     ! energy interval
Nk1 = 61          ! number k points(k点数目增加,计算细节可以体现更好,相应的计算时间也变长)
Nk2 = 101          ! number k points 
NP = 1              ! Number of princple layers for surface green's function(控制表面格林函数计算)
/

LATTICE  ! 元胞信息设置(可参考VASP)
Angstrom
     3.2981915     0.0000000     0.0000000
     0.0000000     5.9443957    -0.0465052
     0.0000000    -0.0987512    28.0482181

ATOM_POSITIONS
6                               ! number of atoms for projectors
Cartisen                          ! Direct or Cartisen coordinate
W        1.6490957     2.3786520    16.4406451
W        0.0000000     0.1332231    16.6311989
Se       1.6490957     0.7737841    18.4325241
Se       0.0000000     3.7846718    17.9156060
Se       0.0000000     1.7381028    14.6392263
Se       1.6490957     4.6716096    15.1096395

PROJECTORS  !设置投影轨道
2*6 4*3        ! number of projectors(每个W原子贡献6个轨道,每个Se原子贡献3个轨道)
W s dz2 dxz dyz dx2-y2 dxy  ! 分别在这里写出是哪些轨道要进行投影(这部分内容我不太懂,还在学习中)
W s dz2 dxz dyz dx2-y2 dxy
Se  pz px py 
Se  pz px py 
Se  pz px py 
Se  pz px py 

SURFACE      !  控制表面态计算
 0  1  0     ! 第一个表面的方向
 0  0  1     ! 第二个表面的方向
 1  0  0     ! 第三个表面的方向

! 当这里计算slab表面态计算的时候,只有前两个矢量方向是有用的,第三个方向默认是开边界的,所以如果向求解不同方向开边界的能谱图只需要对这个参数的顺序进行调整即可,或者就只写出两个周期方向的矢量,比如像看a方向开边界,那么b,c方向就是周期的,这里默认a,b,c是直角坐标的三个基矢.
!--------------------------------------
!SURFACE
!0  1  0
!0  0  1
!---------------------------------------

KPATH_BULK            ! k point path(在进行体态能带计算时候,控制计算的高对称路径)
4              ! number of k line only for bulk band
  X 0.50000  0.00000  0.00000   G   0.00000  0.00000  0.00000   
  G 0.00000  0.00000  0.00000   Y   0.00000  0.50000  0.00000   
  Y 0.00000  0.50000  0.00000   M   0.50000  0.50000  0.00000   
  M 0.50000  0.50000  0.00000   G   0.00000  0.00000  0.00000   

KPATH_SLAB  ! 表面态能谱计算时候的路径选择
2        ! numker of k line for 2D case
X  0.5  0.0 G  0.0  0.0  ! k path for 2D case
G  0.0  0.0 M  0.5  0.5 
! 从上面可以分析出,在进行表面态计算时选择的路径是X--->G--->M

KPLANE_BULK ! 这个参数暂时并不清楚是控制什么计算的
 0.00  0.00  0.00   ! Original point for 3D k plane 
 1.00  0.00  0.00   ! The first vector to define 3d k space plane
 0.00  0.50  0.00   ! The second vector to define 3d k space plane
```

# Results
将上面的两个文件内容准备好之后,就可以进行计算了
```shell
mpirun -np 10 wt.x
```
计算完成之后对结果进行可视化处理
```shell
gnuplot surfdos_l.gnu
```
![png](/assets/images/wannierTools/ti3.png)

# Bi$_2$Se$_3$
这里对3D拓扑绝缘体Bi$_2$Se$_3$进行一些计算

![png](/assets/images/wannierTools/ti4.png)

首先解压紧束缚数据
```shell
tar xzvf wannier90_hr.dat.tar.gz
```
接下来分析一下控制计算的`wt.in`文件
```fortran
&TB_FILE
Hrfile = 'wannier90_hr.dat'      
Package = 'VASP'             ! obtained from VASP, it could be 'VASP', 'QE', 'Wien2k', 'OpenMx'
/

LATTICE
Angstrom
-2.069  -3.583614  0.000000     ! crystal lattice information
 2.069  -3.583614  0.000000
 0.000   2.389075  9.546667

ATOM_POSITIONS
5                               ! number of atoms for projectors
Direct                          ! Direct or Cartisen coordinate
 Bi 0.3990    0.3990    0.6970
 Bi 0.6010    0.6010    0.3030
 Se 0.0000    0.0000    0.5000
 Se 0.2060    0.2060    0.1180
 Se 0.7940    0.7940    0.8820

PROJECTORS
 3 3 3 3 3          ! number of projectors
Bi px py pz         ! projectors
Bi px py pz
Se px py pz
Se px py pz
Se px py pz

SURFACE            ! Specify surface with two vectors, see doc
 1  0  0
 0  1  0


!> bulk band structure calculation flag
&CONTROL
BulkBand_calc            = T ! 计算体态能带
BulkBand_points_calc     = T ! Flag for bulk energy band calculation for some k points
DOS_calc                 = T ! Flag for density of state calculation
SlabBand_calc            = T ! Flag for 2D slab energy band calculation
SlabBandWaveFunc_calc    = T ! Flag for 2D slab band wave function
SlabBand_plane_calc      = T ! Flag for 2D slab energy band calculation
WireBand_calc            = T ! Flag for 1D wire energy band calculation
SlabSS_calc              = T ! Flag for surface state ARPES spectrum calculation
SlabArc_calc             = T ! Flag for surface state fermi-arc calculation
SlabQPI_calc             = T ! Flag for surface state QPI spectrum calculation in a given k plane in 2D BZ
Z2_3D_calc               = T ! Flag for Z2 number calculations of 6 planes(强TI需要有4个indeces来判断)
SlabSpintexture_calc     = T ! Flag for surface state spin-texture calculation(计算表面态上的spin分布变化)
Wanniercenter_calc       = T ! Flag for Wilson loop calculation
/

&SYSTEM
NSLAB = 4               ! Number of slabs for 2d Slab system
NSLAB1= 2               ! Number of slabs for 1D wire system
NSLAB2= 2               ! Number of slabs for 1D wire system 
NumOccupied = 18        !> Number of occupied bands for bulk unit cell
SOC = 1                 ! A parameter to control soc;  Soc=0 means no spin-orbit coupling; Soc>0 means spin-orbit coupling
E_FERMI = 4.4195        ! Fermi energy, search E-fermi in OUTCAR for VASP, set to zero for Wien2k
surf_onsite= 0.0        !> surface onsite energy shift
/

&PARAMETERS
Eta_Arc = 0.001     ! infinite small value, like brodening, used to calculate dos epsilon+i eta
E_arc = 0.0         ! energy level for contour plot of spectrum, Fermi energy for arc calculation
OmegaNum = 400      ! omega number       
OmegaMin = -0.6     ! energy interval
OmegaMax =  0.5     ! energy interval
Nk1 = 101           ! number k points  odd number would be better(如果这个数值过小,结果会非常粗糙)
Nk2 = 101            ! number k points  odd number would be better
Nk3 = 101            ! number k points  odd number would be better
NP = 1              ! number of principle layers
Gap_threshold = 0.01 ! threshold for FindNodes_calc output, threshold value for output the the k points data for Gap3D
/

KPATH_BULK            ! k point path(动量空间能带计算路径)
4              ! number of k line only for bulk band
G 0.00000 0.00000 0.0000 Z 0.00000 0.00000 0.5000
Z 0.00000 0.00000 0.5000 F 0.50000 0.50000 0.0000
F 0.50000 0.50000 0.0000 G 0.00000 0.00000 0.0000
G 0.00000 0.00000 0.0000 L 0.50000 0.00000 0.0000  

KPATH_SLAB   ! 半开边界计算能带时的路径选择
2        ! numker of k line for 2D case
K 0.33 0.67 G 0.0 0.0  ! k path for 2D case
G 0.0 0.0 M 0.5 0.5

KPLANE_SLAB   
-0.1 -0.1      ! Original point for 2D k plane(2D半开边界计算时能带计算的位置选择)
 0.2  0.0      ! The first vector to define 2D k plane 
 0.0  0.2      ! The second vector to define 2D k plane  for arc plots

KPLANE_BULK
 0.00  0.00  0.50   ! Original point for 3D k plane 
 1.00  0.00  0.00   ! The first vector to define 3d k space plane
 0.00  0.50  0.00   ! The second vector to define 3d k space plane


KCUBE_BULK
-0.50 -0.50 -0.50   ! Original point for 3D k plane 
 1.00  0.00  0.00   ! The first vector to define 3d k space plane
 0.00  1.00  0.00   ! The second vector to define 3d k space plane
 0.00  0.00  1.00   ! The third vector to define 3d k cube

WANNIER_CENTRES     ! copy from wannier90.wout
Cartesian
 -0.000040  -1.194745   6.638646 
  0.000038  -1.196699   6.640059 
 -0.000032  -1.192363   6.640243 
 -0.000086  -3.583414   2.908040 
  0.000047  -3.581457   2.906587 
 -0.000033  -3.585864   2.906443 
 -0.000001   1.194527   4.773338 
  0.000003   1.194538   4.773336 
 -0.000037   1.194536   4.773327 
  0.000006  -1.194384   1.130261 
 -0.000018  -1.216986   1.140267 
  0.000007  -1.172216   1.140684 
  0.000011  -3.583770   8.416406 
 -0.000002  -3.561169   8.406398 
 -0.000007  -3.605960   8.405979 
  0.000086  -1.194737   6.638626 
 -0.000047  -1.196693   6.640080 
  0.000033  -1.192286   6.640223 
  0.000040  -3.583406   2.908021 
 -0.000038  -3.581452   2.906608 
  0.000032  -3.585788   2.906424 
  0.000001   1.194548   4.773330 
 -0.000003   1.194537   4.773332 
  0.000037   1.194539   4.773340 
 -0.000011  -1.194381   1.130260 
  0.000002  -1.216981   1.140268 
  0.000007  -1.172191   1.140687 
 -0.000006  -3.583766   8.416405 
  0.000018  -3.561165   8.406400 
 -0.000007  -3.605935   8.405982 
```
将两个主要文件准备好之后,就可以开始计算了
```shell
mpirun -np 2 wt.x &
```
对计算结果进行可视化
```shell
gnuplot bulkek.gnu  # 体态能带计算
gnuplot wanniercenter3D_Z2.gnu-tutorial # Wilson loop计算
gnuplot surfdos_l.gnu  # 表面态密度计算
gnuplot arc_l.gnu  # 计算费米弧
gnuplot spintext_l.gnu  # 计算费米弧上的spin分布
gnuplot slabek.gnu # 计算边界态
```
![png](/assets/images/wannierTools/ti5.png)

![png](/assets/images/wannierTools/ti6.png)

![png](/assets/images/wannierTools/ti7.png)

![png](/assets/images/wannierTools/ti8.png)

上面就是这个实例的一些参数设置和计算得到的一些结果.

# 练习
在这里做个小练习,改动参数来计算一下不同表面上的能带以及谱函数,这里主要通过修改`SURFACE`这个参数下的内容来计算其它表面上的边界态以及能谱,将这个参数修改为
```shell
SURFACE            ! Specify surface with two vectors, see doc
 0  1  0
 0  0  1
 ```
 计算的结果如下

![png](/assets/images/wannierTools/ti9.png)

![png](/assets/images/wannierTools/ti10.png)

![png](/assets/images/wannierTools/ti11.png)

![png](/assets/images/wannierTools/ti12.png)


# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)