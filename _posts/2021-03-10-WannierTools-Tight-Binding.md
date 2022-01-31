---
title: 从紧束缚模型出发构建WannierTools需要的数据
tags: Topology
layout: article
license: true
toc: true
key: a20210310
pageview: true
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
之前学习用WannierTools来计算一些拓扑性质,其中也涉及到了利用它来研究一个紧束缚模型的问题,这里主要是详细分析一下我们如何从紧束缚模型出发来造出WannierTools所需要的数据结构.
<!--more-->
# 紧束缚模型哈密顿量

$$
\begin{equation}
H(\mathbf{k})=A_x \sin(k_x)\sigma_x s_z + A_y \sin(k_y) \sigma_y s_0 + (m_0 - t_x \cos(k_x) - t_y \cos(k_y))\sigma_z s_0 + h_z \sigma_0 s_x\label{ham1} 
\end{equation}
$$

这是一个2D量子自旋霍尔效应, 在其中加上一个Zeeman场$h_z$. 参数选取为$m_0=1.5, t_x=t_y=1.0, A_x=A_y=1.0, h_z=0.2.$ 紧束缚模型哈密顿量的做法就是将三角函数改写成指数函数形式

$$
\begin{equation}
\sin(k_x)=\frac{1}{2i}(e^{ik_x}-e^{-ik_x})\qquad \cos(k_x)=\frac{1}{2}(e^{ik_x}+e^{-ik_x})
\end{equation}
$$

这里假设了晶格常数$a=1$

$$
\begin{equation}
\begin{aligned}
&e^{ik_x}\rightarrow(100)\text{hopping}\qquad e^{-ik_x}\rightarrow(-100)\text{hopping}\\
&e^{ik_y}\rightarrow(010)\text{hopping}\qquad e^{-ik_y}\rightarrow(0-10)\text{hopping}\\
&e^{ik_x+ik_y}\rightarrow(110)\text{hopping}\qquad e^{-ik_x-ik_y}\rightarrow(-1-10)\text{hopping}\\
&e^{ik_x-ik_y}\rightarrow(1-10)\text{hopping}\qquad e^{-ik_x+ik_y}\rightarrow(-110)\text{hopping}\\
\end{aligned}
\end{equation}\label{tri}
$$

接下来可以将(\ref{ham1})中的每一项都改写成关于指数形式的矩阵, 这里矩阵的直积顺序为$s_i\otimes\sigma_i$.

$$
\begin{equation}
{\color{blue}A_x\sin(k_x)\sigma_xs_z=\frac{A_x}{2i}e^{ik_x}
	\left(
	\begin{array}{cccc}
		0 & 1 & 0 & 0 \\
		1 & 0 & 0 & 0 \\
		0 & 0 & 0 & -1 \\
		0 & 0 & -1 & 0 \\
	\end{array}
	\right)-\frac{A_x}{2i}e^{-ik_x}
	\left(
	\begin{array}{cccc}
		0 & 1 & 0 & 0 \\
		1 & 0 & 0 & 0 \\
		0 & 0 & 0 & -1 \\
		0 & 0 & -1 & 0 \\
	\end{array}
	\right)}\label{ax}
\end{equation}
$$

$$
\begin{equation}
{\color{purple}A_y\sin(k_y)\sigma_ys_0=\frac{A_y}{2i}e^{ik_y}\left(
	\begin{array}{cccc}
		0 & -i & 0 & 0 \\
		i & 0 & 0 & 0 \\
		0 & 0 & 0 & -i \\
		0 & 0 & i & 0 \\
	\end{array}
	\right)-\frac{A_y}{2i}e^{-ik_y}\left(
	\begin{array}{cccc}
		0 & -i & 0 & 0 \\
		i & 0 & 0 & 0 \\
		0 & 0 & 0 & -i \\
		0 & 0 & i & 0 \\
	\end{array}
	\right)}\label{ay}
\end{equation}
$$

$$
\begin{equation}
t_x\cos(k_x)\sigma_zs_0=\frac{t_x}{2}e^{ik_x}\left(
\begin{array}{cccc}
	1 & 0 & 0 & 0 \\
	0 & -1 & 0 & 0 \\
	0 & 0 & 1 & 0 \\
	0 & 0 & 0 & -1 \\
\end{array}
\right)+\frac{t_x}{2}e^{-ik_x}\left(
\begin{array}{cccc}
	1 & 0 & 0 & 0 \\
	0 & -1 & 0 & 0 \\
	0 & 0 & 1 & 0 \\
	0 & 0 & 0 & -1 \\
\end{array}
\right)\label{tx}
\end{equation}
$$

$$
\begin{equation}
	t_y\cos(k_y)\sigma_zs_0=\frac{t_y}{2}e^{ik_y}\left(
	\begin{array}{cccc}
		1 & 0 & 0 & 0 \\
		0 & -1 & 0 & 0 \\
		0 & 0 & 1 & 0 \\
		0 & 0 & 0 & -1 \\
	\end{array}
	\right)+\frac{t_y}{2}e^{-ik_y}\left(
	\begin{array}{cccc}
		1 & 0 & 0 & 0 \\
		0 & -1 & 0 & 0 \\
		0 & 0 & 1 & 0 \\
		0 & 0 & 0 & -1 \\
	\end{array}
	\right)\label{ty}
\end{equation}
$$

$$
\begin{equation}
m_0\sigma_zs_0=m_0\left(
\begin{array}{cccc}
	1 & 0 & 0 & 0 \\
	0 & -1 & 0 & 0 \\
	0 & 0 & 1 & 0 \\
	0 & 0 & 0 & -1 \\
\end{array}
\right)\qquad h_z\sigma_0s_x=h_z\left(
\begin{array}{cccc}
	0 & 0 & 1 & 0 \\
	0 & 0 & 0 & 1 \\
	1 & 0 & 0 & 0 \\
	0 & 1 & 0 & 0 \\
\end{array}
\right)\label{on}
\end{equation}
$$

这里的$h_z$与$m_0$中并不涉及hopping, 最后代表的是onsite能量, 即就是$(000)$方向的跃迁.

# 构建WannierTools紧束缚数据
接下来就是结合(\ref{tri},\ref{ax},\ref{ay},\ref{tx},\ref{ty})来构建WannierTools计算所需要的紧束缚模型数据.
其中最主要的结构如下

$$
\begin{equation}
x\quad y\quad z\quad m\quad n\quad\text{Re}\quad\text{Im}
\end{equation}
$$

$$
\begin{equation*}
		\begin{aligned}
			&x\rightarrow\text{$x$方向跃迁取向, 可以取(0,1,-1)}\quad y\rightarrow\text{$y$方向跃迁取向, 可以取(0,1,-1)}\\
			&z\rightarrow\text{$z$方向跃迁取向, 可以取(0,1,-1)}\quad m,n\text{对应的就是具体某一个取向下,对应的矩阵的元素坐索引}\\
			&\text{Re}(m,n)\text{矩阵位置上值的实部}\quad \text{Im}(m,n)\text{矩阵位置上值的实部}
		\end{aligned}
\end{equation*}
$$

下面给一个具体的例子来写这个数据, 比如要写$x$负方向的一个矩阵

$$
\begin{equation}
\frac{A_x}{2i}e^{-ik_x}
\left(
\begin{array}{cccc}
	0 & 1 & 0 & 0 \\
	1 & 0 & 0 & 0 \\
	0 & 0 & 0 & -1 \\
	0 & 0 & -1 & 0 \\
\end{array}
\right)
\end{equation}
$$

$$
\begin{equation}
\left(
\begin{array}{ccccccc}
-1&0&0&1&1&0&0\\
-1&0&0&1&2&0&\frac{A_x}{2i}\\
-1&0&0&1&3&0&0\\
-1&0&0&1&3&0&0\\
-1&0&0&2&1&0&\frac{A_x}{2i}\\
-1&0&0&2&2&0&0\\
-1&0&0&2&3&0&0\\
-1&0&0&2&4&0&0\\
-1&0&0&3&1&0&0\\
-1&0&0&3&2&0&0\\
-1&0&0&3&3&0&0\\
-1&0&0&3&4&0&-\frac{A_x}{2i}\\
-1&0&0&4&1&0&0\\
-1&0&0&4&2&0&0\\
-1&0&0&4&3&0&0\\
-1&0&0&4&4&0&-\frac{A_x}{2i}\\
\end{array}
\right)
\end{equation}
$$

# 产生程序
最后把自己写的构建这个模型数据的Fortran程序附上
```fortran
! Chern insulator
! usage:
! compile and run
! gfortran writeHmnR.f90 -o writehmnr
! ./writehmnr
! H  = Ax sin(kx)sigma_x s_z + Ay sin(ky) sigma_y s_0 + (m0 - tx cos(kx) - ty cos(ky))sigma_z s_0 + lam1 sigma_0 s_x 

program writeHmnR
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
lamc = 0.2d0  ! 层间耦合大小

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

hmnr(1 ,2, ir) = lamc
hmnr(2 ,1, ir) = lamc
hmnr(3 ,4, ir) = lamc
hmnr(4 ,3, ir) = lamc

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
open(unit=105, file='ChernInsulator2L_hr.dat')
write(105, *)'Two Layers Chern Insulator couple'
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

所有这些内容可以[点击这里下载](/assets/pdf/Note-Latex.zip)