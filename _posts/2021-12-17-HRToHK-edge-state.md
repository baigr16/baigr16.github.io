---
title: 利用Wannier90的hr数据计算边界态
tags: Python 
layout: article
license: true
toc: true
key: a20211217
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
这里整理一下通过读取紧束缚模型的hr数据，从而进一步通过之前学习到的迭代格林函数方法来计算边界态，同时也构建一个有限大小的格点来计算边界态，也算是达到对一个模型完整研究的目的。
{:.info}
<!--more-->
# 前言
之前已经学习了如何将一个现有的动量空间中的模型,变换成Wannier90所输出的tight binding的数据,可以参考[从紧束缚模型出发构建WannierTools需要的数据](https://yxli8023.github.io/2021/03/10/WannierTools-Tight-Binding.html)这个Blog中的内容,
因为现有一些程序是可以直接利用这些数据进行计算的,但是同样有时候需要利用实空间的tight binding数据来进行动量空间中的一些计算,这里就想整理一个比较通用的脚本,可以方便的将一个Wannier90输出的TB模型的数据,通过Fourier变换从而来得到动量空间中的哈密顿量,
这样可以方便自己进行后面的一些计算，同时也可以只对某一个方向进行部分Fourier变换，从而构建一个半无限或者有效大小的体系来计算边界态，而半无限大系统计算边界态则是利用迭代格林函数的方法进行计算，有效大小体系的边界态我也会在后面的Blog中整理。
# 数据结构分析
这里就以最熟悉的BHZ模型为例,借用[从紧束缚模型出发构建WannierTools需要的数据](https://yxli8023.github.io/2021/03/10/WannierTools-Tight-Binding.html)这篇Blog中的方法,来构建

$$
H(\mathbf{k})=(m_0-t_x\cos k_x-t_y\cos k_y)\sigma_z+ A_x\sin k_x\sigma_xs_z+A\sin k_y\sigma_y
$$

这个哈密顿量的`wannier90_hr.dat`这个数据，这里说明一下，下面程序中的参数设置可能和上面公式中的有些出入，也就是能带反转是发生在$(\pi,\pi)$点的，所以边界态也同样出现在这个位置，程序如下
```fortran
! Quantum Spin Hall
! usage:
! compile and run
! gfortran hr_bhz.f90 -o hr_bhz
! ./hr_bhz
! H  = Ax sin(kx)sigma_x s_z + Ay sin(ky) sigma_y s_0 + (m0 - tx cos(kx) - ty cos(ky))sigma_z s_0

program hr_bhz
implicit none
integer, parameter :: dp=kind(1d0) ! 双精度计算
complex(dp), parameter :: im = (0d0, 1d0) ! 虚数单位i
complex(dp), parameter :: zzero = (0d0, 0d0)
integer :: i, j
integer :: ir
integer :: nwann
!> arrays for hamiltonian storage
integer :: nrpts
integer, allocatable :: ndegen(:)
integer, allocatable :: irvec(:, :)
complex(dp), allocatable :: hmnr(:, :, :)

!> three lattice constants
real(dp) :: Ax,Ay,m0,tx,ty,lamc
Ax = 1d0
Ay = 1d0
m0 = 1.5d0
tx = 1d0
ty = 1d0

nwann = 2 ! 构建紧束缚模型时候Wannier轨道的数目,这个模型中有4个轨道,化学式为零时占据两个 
nrpts = 17
allocate(irvec(3, nrpts)) ! 该模型时2D的,所以hopping的方向也就只有2个,所以在之后的设置中保持第三个方向上恒为0
allocate(ndegen(nrpts))
allocate(hmnr(nwann*2, nwann*2, nrpts))
irvec = 0
ndegen = 1 ! 手动设置紧束缚模型数据所需要
hmnr = zzero ! 矩阵初始化


! 0 0 onsite能量
ir = 1  ! ir用来记录这里一共会有多少个hopping的方位(包括了onsite)
irvec(1, ir) =  0
irvec(2, ir) =  0
hmnr(1, 1, ir) = m0 
hmnr(2, 2, ir) = m0
hmnr(3, 3, ir) = -m0
hmnr(4, 4, ir) = -m0


!1 0 x正反向hopping
ir = ir + 1
irvec(1, ir) =  1
irvec(2, ir) =  0


hmnr(1, 1, ir) = tx/2.0
hmnr(2, 2, ir) = tx/2.0
hmnr(3, 3, ir) = -tx/2.0
hmnr(4, 4, ir) = -tx/2.0

hmnr(1, 3, ir) = ax/(2.0*im)
hmnr(2, 4, ir) = -ax/(2.0*im)
hmnr(3, 1, ir) = ax/(2.0*im)
hmnr(4, 2, ir) = -ax/(2.0*im)


!-1 0 x负方向hopping
ir = ir+ 1
irvec(1, ir) = -1
irvec(2, ir) =  0

hmnr(1, 1, ir) = tx/2.0
hmnr(2, 2, ir) = tx/2.0
hmnr(3, 3, ir) = -tx/2.0
hmnr(4, 4, ir) = -tx/2.0

hmnr(1, 3, ir) = -ax/(2.0*im)
hmnr(2, 4, ir) = ax/(2.0*im)
hmnr(3, 1, ir) = -ax/(2.0*im)
hmnr(4, 2, ir) = ax/(2.0*im)

! 0 1 y正方向hopping
ir = ir + 1
irvec(1, ir) =  0 
irvec(2, ir) =  1

hmnr(1, 1, ir) = ty/2.0
hmnr(2, 2, ir) = ty/2.0
hmnr(3, 3, ir) = -ty/2.0
hmnr(4, 4, ir) = -ty/2.0

hmnr(1, 3, ir) = -im*ay/(2.0*im)
hmnr(2, 4, ir) = -im*ay/(2.0*im)
hmnr(3, 1, ir) = im*ay/(2.0*im)
hmnr(4, 2, ir) = im*ay/(2.0*im)

!0 -1  y负方向hopping
ir = ir + 1
irvec(1, ir) =  0 
irvec(2, ir) = -1

hmnr(1, 1, ir) = ty/2.0
hmnr(2, 2, ir) = ty/2.0
hmnr(3, 3, ir) = -ty/2.0
hmnr(4, 4, ir) = -ty/2.0

hmnr(1, 3, ir) = im*ay/(2.0*im)
hmnr(2, 4, ir) = im*ay/(2.0*im)
hmnr(3, 1, ir) = -im*ay/(2.0*im)
hmnr(4, 2, ir) = -im*ay/(2.0*im)

nrpts = ir

!> write to new_hr.dat
open(unit=105, file='wannier90_hr.dat')
write(105, *)'BHZ model'
write(105, *)nwann*2
write(105, *)nrpts
write(105, '(15I5)')(ndegen(i), i=1, nrpts)
do ir = 1, nrpts
   do i = 1, nwann*2
      do j = 1, nwann*2
         write(105, '(5I5, 2f16.8)')irvec(:, ir), i, j, HmnR(i, j, ir)
      end do
   end do
end do
close(105)
stop
end ! end of program
```
先来看看数据结构

![png](/assets/images/wannier90/hr-1.png)

- 第一行
> 注释,说明一下这个TB数据是什么体系的
- 第二行
> 能带数目,从模型的角度可以理解成哈密顿量的维度,多大的矩阵自然就对应几条能带,但是可能包含简并,简并度是多少就算多少条能带.
- 第三行
> 在实空间中,每个格点上的hopping的方向,比如$\cos k_x=\frac{1}{2}(e^{ik_x\cdot R_x}+e^{-ik_x\cdot R_x})\quad R_x=1$,这里假设晶格长度为单位1,所以此时就是有两个hopping的方向,同样的格点上的onsite占位能也同样算是一个方向,比如在BHZ模型中$m_0$就代表的是格点上的占位能,它同样也是一个方向,所以BHZ中x的正负方向,以及y的正负方向,还有onsite能量,就一共存在5个hopping的方向,这也就是第三行的含义.
- 第四行
> 这一行是每个hopping方向格点上的简并度,因为这里研究的是BHZ模型,TB数据是手动构造的,所以这个数字都是1,这里由多少个hopping方向就会存在所少个数字.
- 第五行及之后

$$
x\quad y\quad z\quad m\quad n\quad\text{Re}\quad\text{Im}
$$

这里前面的三个数(x,y,z)就是用来定义实空间中hopping的方向,比如x正方向最近邻hopping为$(1,0,0)$,其他更远或者其他方向的hopping就可以通过类似的方式来得到.后面的$(m,n)$表示的是能带指标,而最后两行代表的是hopping对应的实部和虚部.关于数据结构这部分的内容,同样可以参考[WannierTools](http://www.wanniertools.com/input.html)中的解释.
```fortran
 BHZ model
           4
           5
    1    1    1    1    1
    0    0    0    1    1      1.50000000      0.00000000
    0    0    0    1    2      0.20000000      0.00000000
    0    0    0    1    3      0.00000000      0.00000000
    0    0    0    1    4      0.00000000      0.00000000
    0    0    0    2    1      0.20000000      0.00000000
    0    0    0    2    2      1.50000000      0.00000000
    0    0    0    2    3      0.00000000      0.00000000
    0    0    0    2    4      0.00000000      0.00000000
    0    0    0    3    1      0.00000000      0.00000000
    0    0    0    3    2      0.00000000      0.00000000
    0    0    0    3    3     -1.50000000      0.00000000
    0    0    0    3    4      0.20000000      0.00000000
    0    0    0    4    1      0.00000000      0.00000000
    0    0    0    4    2      0.00000000      0.00000000
    0    0    0    4    3      0.20000000      0.00000000
    0    0    0    4    4     -1.50000000      0.00000000
    1    0    0    1    1      0.50000000      0.00000000
    1    0    0    1    2      0.00000000      0.00000000
    1    0    0    1    3      0.00000000     -0.50000000
    1    0    0    1    4      0.00000000      0.00000000
    1    0    0    2    1      0.00000000      0.00000000
    1    0    0    2    2      0.50000000      0.00000000
    1    0    0    2    3      0.00000000      0.00000000
    1    0    0    2    4      0.00000000      0.50000000
    1    0    0    3    1      0.00000000     -0.50000000
    1    0    0    3    2      0.00000000      0.00000000
    1    0    0    3    3     -0.50000000      0.00000000
    1    0    0    3    4      0.00000000      0.00000000
    1    0    0    4    1      0.00000000      0.00000000
    1    0    0    4    2      0.00000000      0.50000000
    1    0    0    4    3      0.00000000      0.00000000
    1    0    0    4    4     -0.50000000      0.00000000
   -1    0    0    1    1      0.50000000      0.00000000
   -1    0    0    1    2      0.00000000      0.00000000
   -1    0    0    1    3      0.00000000      0.50000000
   -1    0    0    1    4      0.00000000      0.00000000
   -1    0    0    2    1      0.00000000      0.00000000
   -1    0    0    2    2      0.50000000      0.00000000
   -1    0    0    2    3      0.00000000      0.00000000
   -1    0    0    2    4      0.00000000     -0.50000000
   -1    0    0    3    1      0.00000000      0.50000000
   -1    0    0    3    2      0.00000000      0.00000000
   -1    0    0    3    3     -0.50000000      0.00000000
   -1    0    0    3    4      0.00000000      0.00000000
   -1    0    0    4    1      0.00000000      0.00000000
   -1    0    0    4    2      0.00000000     -0.50000000
   -1    0    0    4    3      0.00000000      0.00000000
   -1    0    0    4    4     -0.50000000      0.00000000
    0    1    0    1    1      0.50000000      0.00000000
    0    1    0    1    2      0.00000000      0.00000000
    0    1    0    1    3     -0.50000000      0.00000000
    0    1    0    1    4      0.00000000      0.00000000
    0    1    0    2    1      0.00000000      0.00000000
    0    1    0    2    2      0.50000000      0.00000000
    0    1    0    2    3      0.00000000      0.00000000
    0    1    0    2    4     -0.50000000      0.00000000
    0    1    0    3    1      0.50000000      0.00000000
    0    1    0    3    2      0.00000000      0.00000000
    0    1    0    3    3     -0.50000000      0.00000000
    0    1    0    3    4      0.00000000      0.00000000
    0    1    0    4    1      0.00000000      0.00000000
    0    1    0    4    2      0.50000000      0.00000000
    0    1    0    4    3      0.00000000      0.00000000
    0    1    0    4    4     -0.50000000      0.00000000
    0   -1    0    1    1      0.50000000      0.00000000
    0   -1    0    1    2      0.00000000      0.00000000
    0   -1    0    1    3      0.50000000      0.00000000
    0   -1    0    1    4      0.00000000      0.00000000
    0   -1    0    2    1      0.00000000      0.00000000
    0   -1    0    2    2      0.50000000      0.00000000
    0   -1    0    2    3      0.00000000      0.00000000
    0   -1    0    2    4      0.50000000      0.00000000
    0   -1    0    3    1     -0.50000000      0.00000000
    0   -1    0    3    2      0.00000000      0.00000000
    0   -1    0    3    3     -0.50000000      0.00000000
    0   -1    0    3    4      0.00000000      0.00000000
    0   -1    0    4    1      0.00000000      0.00000000
    0   -1    0    4    2     -0.50000000      0.00000000
    0   -1    0    4    3      0.00000000      0.00000000
    0   -1    0    4    4     -0.50000000      0.00000000
```

# 文件读写
```python
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------
file = open('wannier90_hr.dat') # 打开TB模型的数据
comment = file.readline() 
Band_number = int(file.readline())
Hopping_number = int(file.readline()) # 对于单独的一个字符串可以将其转化为整数

degenerate = [] # 用来存储每个hopping的简并度,每个hopping方向对应一个数字
con = 0
for line in file: # 接着上面读取到行的位置继续读
    temp = line.strip().split() # 移除字符串头尾指定的字符（默认为空格）,并以空格将这些字符串分开
    t1 = [int(x) for x in temp] #遍历每个字符串并转换为整数
    degenerate.extend(t1) # 数据追加
    if len(degenerate) == Hopping_number: # 当获取的hopping对应的简并度的个数和hopping方向的数目相等的时候,就停止读取文件
        break
# 接下来读取核心的hr参数,也就是不同带之间hopping的大小以及方向,hr文件中剩下所有的数据都是这些信息,因此可以直接读取到文件末尾
Hr = []
for line in file:
    temp = line.strip().split() # 先除去空格,再以空格分立成单独的字符串
    # 将所有的字符串再转化成数值进行存储
    Hr.append([int(temp[0]),int(temp[1]),int(temp[2]),int(temp[3]),int(temp[4]),float(temp[5]),float(temp[6])])
file.close() # 读取完成之后关闭文件
Hr1 = np.array(Hr) # 将一个list转化成矩阵

# 在通过TB数据进行Fourier变换得到动量空间中的模型时,需要利用的hopping R,这里将这些hopping方向R单独提取出来
Hopping_dir = np.zeros([Hopping_number,3],np.int64)
for i in range(Hopping_number):
    Hopping_dir[i] = Hr1[Band_number**2 * i ,0:3] # 获取所有的hopping方向
```
## 规范一

$$
\rvert\chi_i^\mathbf{k}\rangle=\sum_\mathbf{R}e^{i\mathbf{k}\cdot(\mathbf{R}-\mathbf{\tau}_i)}\rvert\phi_{\mathbf{R}i}\rangle
$$

利用其构建Bloch态

$$
\rvert\psi_{n\mathbf{k}}\rangle=\sum_jC_{j}^{n\mathbf{k}}\rvert\chi_j^\mathbf{k}\rangle
$$

$$
H_{ij}^\mathbf{k}=\langle\chi_i^\mathbf{k}\rvert H\rvert\chi_j^\mathbf{k}\rangle=\sum_\mathbf{R}e^{i\mathbf{k}\cdot(\mathbf{R}+\tau_j-\tau_i)}H_{ij}(\mathbf{R})
$$

$$
H(\mathbf{k})C_{n\mathbf{k}}=E_{n\mathbf{k}}C_{n\mathbf{k}}
$$

这里的$ij$是轨道(能带)的索引.

## 规范二

$$
\rvert\tilde{\chi}_i^\mathbf{k}\rangle=\sum_\mathbf{R}e^{i\mathbf{k}\cdot\mathbf{R}}\rvert\phi_{\mathbf{R}i}\rangle
$$

利用其构建Bloch态

$$
\rvert\tilde{\psi}_{n\mathbf{k}}\rangle=\sum_j\tilde{C}_{j}^{n\mathbf{k}}\rvert\tilde{\chi}_j^\mathbf{k}\rangle
$$

$$
\tilde{H}_{ij}^\mathbf{k}=\langle\tilde{\chi}_i^\mathbf{k}\rvert H\rvert\tilde{\chi}_j^\mathbf{k}\rangle=\sum_\mathbf{R}e^{i\mathbf{k}\cdot\mathbf{R}}H_{ij}(\mathbf{R})
$$

其满足的本征方程为

$$
\tilde{H}(\mathbf{k})\tilde{C}_{n\mathbf{k}}=E_{n\mathbf{k}}\tilde{C}_{n\mathbf{k}}
$$

在同TB构建Bloch哈密顿量的时候会有两种规范选择,其实不同之处就在于构建的时候,在hopping中是否考虑轨道的中心$\tau_i$,而通过这两种规范得到的Bloch哈密顿量之间的关系为

$$
\tilde{H}(\mathbf{k})_{ij}=e^{i\mathbf{k}\cdot(\tau_i-\tau_j)}H_{ij}(\mathbf{k})
$$

其实这个额外的位相因子通常是反映在本征矢量上的

$$
\tilde{C}_j^{n\mathbf{k}}=e^{i\mathbf{k}\cdot\tau_j}C^{n\mathbf{k}}_j
$$

这里的n代表的是哈密顿量的那个本征值,$j$则是轨道(能带)的索引.

因为我在这里考虑的是BHZ模型,所以就直接采用第二种规范,不考虑每个轨道的$\tau_i$,主要是对于一个给出形式的Bloch哈密顿量,也确实不知道它每个给定的轨道对应的$\tau_i$是多少,但是如果是利用Wannier90计算出的hr数据(具体材料),那么在hr中就会有对应的每个轨道的$\tau_i$,这个时候就可以同时用两种不同的规范来进行Bloch哈密顿量构建.
```python
#-------------------------------------------------------------------
def Hamk(hamhr,Hopping_dir,Hopping_number,Band_number,kx,ky):
    # 同过TB数据构建动量空间哈密顿量
    H = np.zeros([Band_number,Band_number],np.complex128) # 动量空间哈密顿量
    for i0 in range(Hopping_number):# 对所有的R求和进行Fourier变换
        H += np.exp(-1j*np.dot(Hopping_dir[i0,0:2],[kx,ky]))*hamhr[i0,:,:] # 每个轨道index对应的前面的位相因子是相同的
    return H
```
从上面的公式中可以看到这里需要的就只是不同轨道（能带）之间的hopping大小$H_{ij}$以及hopping的方向$\mathbf{R}$，我们这里采用的是第二种规范，所以实际上就是对所有的$\mathbf{R}$求和，这里每个轨道自由度对应的这个位相因子是相同的，所以在程序中直接就使用了**hamhr[i0,:,:]**，整个这部分用到的核心公式就是

$$
H_{mn}(\mathbf{k})=\sum_\mathbf{R}e^{i\mathbf{k}\cdot\mathbf{R}}H_{mn}(\mathbf{R})
$$

通过上面的方式就完全由TB数据构建出了动量空间中的哈密顿量。这里想要计算边界态，那么在进行Fourier变换的时候就不是要对所有的方向都进行，而只需要对单独的一个方向进行Fourier变换就可以，我在这里想要令$k_y$是个好量子数，而$x$方向仍然是在实空间中的，那么上面的变换公式就需要稍微改动一下

$$H_{mn}(k_y)=\sum_{R_y}e^{ik_yR_y}H_{mn}(R_y)$$

其实这里本质上也就是在之前的基础上，只变换一个方向而已，那么这个变换在程序上的实现如下

```python
def Hamk(hamhr,Hopping_dir,Hopping_number,Band_number,ky):
    # 同过TB数据构建动量空间哈密顿量，这里沿着y方向变换到动量空间，x方向仍然在实空间
    H = np.zeros([Band_number,Band_number],np.complex128) # 动量空间哈密顿量
    for i0 in range(Hopping_number): # 对所有的R求和进行Fourier变换        
        if Hopping_dir[i0,0] != 0& Hopping_dir[i0,1] == 0:
            H += 0 # 这里在碰到x方向的hopping的时候，从公式变换中将这种情况去掉，它是不贡献的，只有ky方向进行Fourier变换
        else:
            H += np.exp(-1j*np.dot(Hopping_dir[i0,0:2],[0,ky]))*hamhr[i0,:,:] # 每个轨道index对应的前面的位相因子是相同的
    return H,hamhr[1] # 返回onsite项与hopping项
```
通过上面的这个函数可以得到元胞内的hopping`H00`和元胞之间的hopping`H01`,那么就可以沿着开边界的方向构建一个有限大小的格点来计算边界态
```python
def cylinder(ky):
    # 利用tb数据构造一个有限大小的体系
    on = 10 # 实空间格点数目
    hamhr,Hopping_dir,Band_number,Hopping_number = Hamhr()# 读取hr中的数据，获得实空间hopping信息
    H00,H01 = Hamk(hamhr,Hopping_dir,Hopping_number,Band_number,ky)
    ham = np.zeros((on*Band_number,on*Band_number),dtype = complex)
    
    
    for i0 in range(0,on): # index for lattice 
        if i0 == 0: # 右边界只有向右hopping
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
        elif i0 == on - 1:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
        else:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
    return ham
```
通过对角化这个哈密顿量就可以得到边界态.

同样的,同样可以通过上面的方式得到构建好的lattice之后,从其中单独再获取元胞内的hopping`H00`和元胞之间的hopping`H01`
```python
def cylinderToGF(ky):
    # 利用tb数据构造一个有限大小的体系
    on = 10 # 实空间格点数目
    hamhr,Hopping_dir,Band_number,Hopping_number = Hamhr()# 读取hr中的数据，获得实空间hopping信息
    H00,H01 = Hamk(hamhr,Hopping_dir,Hopping_number,Band_number,ky)
    ham = np.zeros((on*Band_number,on*Band_number),dtype = complex)
    for i0 in range(0,on): # index for lattice 
        if i0 == 0: # 右边界只有向右hopping
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
        elif i0 == on - 1:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
        else:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
    #-----------------------------------------------------------
    # 在构建好了cylinder之后,就可以直接得到H00和H01,这里使用的是BHZ模型,所以只是包含了最近邻的hopping,
    H00 = np.zeros((Band_number,Band_number),dtype = complex)
    H01 = np.zeros((Band_number,Band_number),dtype = complex)
    for i0 in range(Band_number):
        for i1 in range(Band_number):
            H00[i0,i1] = ham[i0,i1]
            H01[i0,i1] = ham[i0,Band_number + i1]
    return H00,H01
```
这个时候就又可以通过迭代格林函数计算表面态了.
# 代码
这里整理一下完整的代码
```python
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import seaborn as sns
#------------------------------------------------------------------------------------
def Hamhr():
    # ------------------------------------
    file = open('wannier90_hr.dat') # 打开TB模型的数据(通常默认文件名为wannier90_hr.dat)
    comment = file.readline() 
    Band_number = int(file.readline())
    Hopping_number = int(file.readline()) # 对于单独的一个字符串可以将其转化为整数

    degenerate = [] # 用来存储每个hopping的简并度,每个hopping方向对应一个数字
    con = 0
    for line in file: # 接着上面读取到行的位置继续读
        temp = line.strip().split() # 移除字符串头尾指定的字符（默认为空格）,并以空格将这些字符串分开
        t1 = [int(x) for x in temp] #遍历每个字符串并转换为整数
        degenerate.extend(t1) # 数据追加
        if len(degenerate) == Hopping_number: # 当获取的hopping对应的简并度的个数和hopping方向的数目相等的时候,就停止读取文件
            break
    # 接下来读取核心的hr参数,也就是不同带之间hopping的大小以及方向,hr文件中剩下所有的数据都是这些信息,因此可以直接读取到文件末尾
    Hr = []
    for line in file:
        temp = line.strip().split() # 先除去空格,再以空格分立成单独的字符串
        # 将所有的字符串再转化成数值进行存储
        Hr.append([int(temp[0]),int(temp[1]),int(temp[2]),int(temp[3]),int(temp[4]),float(temp[5]),float(temp[6])])
    file.close() # 读取完成之后关闭文件
    Hr = np.array(Hr) # 将一个list转化成矩阵, 这里存储的是核心的hooping方向以及大小

    # 在通过TB数据进行Fourier变换得到动量空间中的模型时,需要利用的hopping R,这里将这些hopping方向R单独提取出来
    Hopping_dir = np.zeros([Hopping_number,3],np.int64)
    for i in range(Hopping_number):
        Hopping_dir[i] = Hr[Band_number**2 * i ,0:3] # 获取所有的hopping方向

    #现在得到了每个轨道hopping方向和对应的大小,那么将其全部存储,方便下一步构建Bloch哈密顿量

    con1 = 0 # 用来计数,工具人
    ham = np.zeros([Hopping_number, Band_number, Band_number], np.complex128)
    for ih in range(Hopping_number):
        for i1 in range(Band_number):
            for i2 in range(Band_number):
                ham[ih,i1,i2] = Hr[con1,5] + 1j*Hr[con1,6] # 对每一个hopping方向,存储每个轨道之间的hopping的大小
                con1 = con1 + 1
    return ham,Hopping_dir,Band_number,Hopping_number
#-----------------------------------------------------------------------------------------
def Hamk(hamhr,Hopping_dir,Hopping_number,Band_number,ky):
    # 同过TB数据构建动量空间哈密顿量，这里沿着y方向变换到动量空间，x方向仍然在实空间
    H = np.zeros([Band_number,Band_number],np.complex128) # 动量空间哈密顿量
    for i0 in range(Hopping_number): # 对所有的R求和进行Fourier变换        
        if Hopping_dir[i0,0] != 0& Hopping_dir[i0,1] == 0:
            H += 0 # 这里在碰到x方向的hopping的时候，从公式变换中将这种情况去掉，它是不贡献的，只有ky方向进行Fourier变换
        else:
            H += np.exp(-1j*np.dot(Hopping_dir[i0,0:2],[0,ky]))*hamhr[i0,:,:] # 每个轨道index对应的前面的位相因子是相同的
    return H,hamhr[1] # 返回onsite项与hopping项
#----------------------------------------------------------------------------------------
def cylinder(ky):
    # 利用tb数据构造一个有限大小的体系
    on = 10 # 实空间格点数目
    hamhr,Hopping_dir,Band_number,Hopping_number = Hamhr()# 读取hr中的数据，获得实空间hopping信息
    H00,H01 = Hamk(hamhr,Hopping_dir,Hopping_number,Band_number,ky)
    ham = np.zeros((on*Band_number,on*Band_number),dtype = complex)
    
    
    for i0 in range(0,on): # index for lattice 
        if i0 == 0: # 右边界只有向右hopping
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
        elif i0 == on - 1:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
        else:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
    return ham
#-------------------------------------------------------------------------------------------------------------
def cylinderToGF(ky):
    # 利用tb数据构造一个有限大小的体系
    on = 10 # 实空间格点数目
    hamhr,Hopping_dir,Band_number,Hopping_number = Hamhr()# 读取hr中的数据，获得实空间hopping信息
    H00,H01 = Hamk(hamhr,Hopping_dir,Hopping_number,Band_number,ky)
    ham = np.zeros((on*Band_number,on*Band_number),dtype = complex)
    for i0 in range(0,on): # index for lattice 
        if i0 == 0: # 右边界只有向右hopping
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
        elif i0 == on - 1:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
        else:
            for i1 in range(Band_number):
                for i2 in range(Band_number):
                    ham[i0*Band_number + i1,i0*Band_number + i2] = H00[i1,i2]  # 好量子数方向
                    ham[i0*Band_number + i1,(i0 + 1)*Band_number + i2] = H01[i1,i2]# 实空间hopping
                    ham[i0*Band_number + i1,(i0 - 1)*Band_number + i2] = np.conj(H01[i1,i2])# 实空间hopping
    #-----------------------------------------------------------
    # 在构建好了cylinder之后,就可以直接得到H00和H01,这里使用的是BHZ模型,所以只是包含了最近邻的hopping,
    H00 = np.zeros((Band_number,Band_number),dtype = complex)
    H01 = np.zeros((Band_number,Band_number),dtype = complex)
    for i0 in range(Band_number):
        for i1 in range(Band_number):
            H00[i0,i1] = ham[i0,i1]
            H01[i0,i1] = ham[i0,Band_number + i1]
    return H00,H01
#-----------------------------------------------------------------------------------------
def Iteration(omega,ki):
    err = 1e-16
    eta = 0.01
    iternum = 200
#     H00 = np.zeros((2,2),np.complex128)
#     H01 = np.zeros((2,2),np.complex128)
#     H00,H01 = hamset(ki)
    H00,H01 = cylinderToGF(ki)
    epsiloni = H00
    epsilons = H00
    epsilons_t = H00
    alphai = H01
    betai = H01.T.conjugate() # 转置共轭
    omegac = omega + eta*1j
    hn = H00.shape[0]
    s0 = np.eye(hn,dtype = complex)
    for i0 in range(iternum):
        g0dem = omegac*s0 - epsiloni
        
        g0 = np.linalg.inv(g0dem)
        
        mat1 = np.dot(alphai,g0)
        
        mat2 = np.dot(betai,g0)

        g0 = np.dot(mat1,betai)
        
        epsiloni = epsiloni + g0

        epsilons = epsilons + g0

        g0 = np.dot(mat2,alphai)

        epsiloni = epsiloni + g0

        epsilons_t = epsilons_t + g0

        g0 = np.dot(mat1, alphai)
        
        alphai = g0

        g0 = np.dot(mat2,betai)

        betai = g0
          
        real_temp = np.sum(np.concatenate(np.abs(alphai)))
        
        if (real_temp < err):
            break
        
    GLLdem = omegac*s0 - epsilons
    GLL = np.linalg.inv(GLLdem)
    #  GLL = epsilons
    GLL =  np.sum(np.concatenate(np.abs(GLL)))


    GRRdem = omegac*s0 - epsilons_t
    GRR = np.linalg.inv(GRRdem)
    GRR =  np.sum(np.concatenate(np.abs(GRR)))
    #  GRR = epsilons_t

    GBdem = omegac*s0 - epsiloni
    GB = np.linalg.inv(GBdem)
    GB =  np.sum(np.concatenate(np.abs(GB)))
    #  GB = epsiloni
    
    return GLL,GRR,GB
#------------------------------------------------------------
def surface():
    nx = 100
    max_omg = 1.5
    re = np.zeros((len(range(-nx,nx))**2,5))
    con = 0
    ix = -1
    iy = -1
    GLL = np.zeros((len(range(-nx,nx)),len(range(-nx,nx))))
    GRR = np.zeros((len(range(-nx,nx)),len(range(-nx,nx))))
    GB = np.zeros((len(range(-nx,nx)),len(range(-nx,nx))))
    for i0 in range(-nx,nx):
        kx = np.pi*i0/nx
        for i1 in range(-nx,nx):
            omg = max_omg*i1/nx
            re1,re2,re3 = Iteration(omg,kx)
            re[con,0] = kx
            re[con,1] = omg
            re[con,2] = re1
            re[con,3] = re2
            re[con,4] = re3
            
            GLL[iy,ix] = np.log(re1)
            GRR[iy,ix] = np.log(re2)
            GB[iy,ix] = np.log(re3)
            con += 1
            iy += 1
        ix += 1
        iy = 0
#     np.savetxt("GLL.dat", [kilist,re1list], fmt="%15.10f")
#     np.savetxt("GRR.dat", [kilist,re1list], fmt="%15.10f")
#     np.savetxt("GB.dat", [kilist,re1list], fmt="%15.10f")
    np.savetxt("density.dat",re , fmt="%15.10f")
    return GLL,GRR,GB
#-------------------------------------------------------------------------------------------------------------
def ishermitian(ham):
    xn,yn = ham.shape 
    sign = 0
    for i0 in range(xn):
        for i1 in range(yn):
            if ham[i0,i1] != np.conj(ham[i1,i0]):
                sign = 1
                print("Hamiltonian isn't Hermitian")
            else:
                sign = 0
    return sign
#-----------------------------------------------------------------------------------
def PltEdgeState():
    # 利用有限大小体系绘制边界态
    print("Start with finite Lattice mode")
    tstart = time.time()
    nky = 100
    hn = (cylinder(0.0)).shape[0]
    kylist = np.linspace(-np.pi, np.pi, nky)
    vallist = []
    for ky in kylist:
        ham = cylinder(ky)
        ishermitian(ham)
        val,vec = np.linalg.eigh(ham)
        vallist.append(np.real(val))
    plt.plot(kylist,vallist)
    plt.show()
    tend = time.time()
    print("Lattice method use time is %.5f" % (tend - tstart))
    # 利用迭代格林函数绘制边界态
    print("Start with Iteration method")
    tstart = time.time()
    GLL,GRR,GB = surface()
    tend = time.time()
    # 绘图
    sns.set()
    ax = sns.heatmap(GRR)
    plt.show()
    print("Iteration time is %.5f" % (tend - tstart))
# -----------------------------------------------------------------------------------
def main():
    PltEdgeState()
#------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
    
```
这里就同时用两种方法计算了边界态.
![png](/assets/images/wannier90/hr-4.png)


# 参考
- 1. Berry Phases in Electronic Structure Theory Electric Polarization, Orbital Magnetization and Topological Insulators

