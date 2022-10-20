---
title: 边界极化
tags:  Latex
layout: article
license: true
toc: true
key: a20221011
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
在计算发生边界拓扑相变的高阶拓扑绝缘体的时候, 有时候会遇到需要计算体系的边界极化, 其实也就是要在一个cylinder结构上计算其对应的Wilson loop, 并计算每个格点上对极化的贡献, 这里就先详细的整理一下理论背景, 后面给一个详细的例子来给出计算记过.
{:.info}
<!--more-->
# 边界极化
首先考虑一个2D系统$\mathbf{k}=(k_x,k_y)$, 将其一个方向变换到实空间

$$c_{k,\alpha}\rightarrow c_{k_x,R_y,\alpha}$$

此时二次量子化形式下的哈密顿量为

$$\begin{align}
H = \sum_{k_x} c^\dagger_{k_x,R_y,\alpha} [h_{k_x}]^{R_y,\alpha,R'_y,\beta} c_{k_x,R'_y,\beta},
\label{eq1}
\end{align}$$

这里的$\alpha,\beta\in 1,\cdots N_{\rm orb}$为轨道指标, $R_y.R_y'\in 1\cdots N$则是$y$方向开边界之后的格点指标, 因为此时$x$方向仍然是周期的, 所以$k_x$还是一个好量子数, 此时计算需要的哈密顿量矩阵为

$$\begin{align}
[h_{k_x}]^{R_y,\alpha,R'_y,\beta} = \sum_n [u^n_{k_x}]^{R_y,\alpha} \epsilon_{n,k_x} [u^{*n}_{k_x}]^{R'_y\beta},
\label{eq2}
\end{align}$$

这里的求和指标$n\in 1\cdots N_{\rm orb}\times N_y$. 如果二维系统的$x$方向和$y$方向都是周期的, 那么$h(k_x,k_y)$是有$N_{\rm orb}$个轨道的哈密顿量, 那么此时方程\eqref{eq2}中$h(k_x)$描述的就是一个具有$N_{\rm orb}\times R_y$个占据态能带的准1D体系, 与[重学Wilson loop ](https://yxli8023.github.io/2022/10/10/Wilsonloop-restudy.html)中一样, 此时将哈密顿量进行对角化

$$
\begin{align}
H = \sum_{n,{k_x}} \gamma^\dagger_{n,{k_x}} \epsilon_{n,k_x} \gamma_{n,{k_x}},
\end{align}
$$

此时准粒子算符与电子算符之间的变换关系为

$$\begin{align}
\gamma_{n,k_x} = \sum_{R_y,\alpha} [u^{*n}_{k_x}]^{R_y,\alpha} c_{k_x,R_y,\alpha}.
\label{eq4}
\end{align}$$

这里就直接利用[重学Wilson loop ](https://yxli8023.github.io/2022/10/10/Wilsonloop-restudy.html)中计算Wilson loop的方法了, 可以得到

$$\begin{align}
[G_{k_x}]^{mn} \equiv \sum_{R_y,\alpha}[u^{*m}_{k_x+\Delta k_x}]^{R_y,\alpha} [u^n_{k_x}]^{R_y,\alpha},
\end{align}$$

此时虽然沿着$y$方向是开放边界, 但是$k_x$还是好量子数, 所以还是可以沿着$k_x$方向计算Wilson loop的, 只不过此时占据态的数量为为$N_y\times N_{\rm orb}$, 而且对应的可以得到此时的Hybird Wannier function

$$\begin{align}
\ket{\Psi^j_{R_x}} = \frac{1}{\sqrt{N_x}} \sum_{n=1}^{N_{occ} \times N_y}\sum_{k_x} \left[ \nu^j_{k_x} \right]^n e^{-i k_x R_x} \gamma^\dagger_{n,k_x}\ket{0},
\label{eq3}
\end{align}$$

这里$j\in 1\cdots N_{\rm occ}\times N_y, R_x\in 1\cdots N_x$, $[v_{k_x}^j]^n$是第$j$个Wilson loop本征态$\rvert v_{k_x}^j\rangle$的第$n$个分量, $\gamma^\dagger_{n,k_x}$如公式\eqref{eq4}所示. 在得到了Hybrid Wannier函数之后, 就可以得到其对应的几率密度

$$\begin{align}
\rho^{j,R_x}(R_y) &= \sum_{R'_x,\alpha} \braket{\Psi^j_{R_x}}{\phi^{R_y, \alpha}_{R'_x}}\braket{\phi^{R_y, \alpha}_{R'_x}}{\Psi^j_{R_x}}\nonumber\\
&=\frac{1}{N_x} \sum_{k_x, \alpha} \left| [u^n_{k_x}]^{R_y, \alpha}[\nu^j_{k_x}]^n\right|^2
\end{align}$$

通过$\rho^{j,R_x}(R_y)$就可以进一步在$y$方向的格点上解析出极化在$x$方向的分量.

特别有用的一点就是可以得到在边界上$R_y=1,N_y$,是否存在具有的Hybrid Wannier函数. 通过几率密度$\rho^{j,R_x}(R_y)$可以计算极化在$x$方向的分量

$$\begin{align}
p_x(R_y) = \sum_j \rho^j(R_y) \nu_x^j 
\label{eq5}
\end{align}$$

从而就可以得到极化在$y$边界上是如何分布的, 通过变化格点$R_y$即可. 为了得到边界极化, 可以从边界中点到边界的上对极化分布进行求和

$$p_x^{\rm edge, y}=\sum_{i=1}^{N_y/2}p_x(i_y)$$

从而就可以得到边界极化. 上面的这些过程同样也可以是沿着$y$方向开边界进行计算的

# 参考
- [Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)
