---
title: 利用投影子方法求解Square lattice边界态
tags: Topology Math Study
layout: article
license: true
toc: true
key: a20211121a
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
在Square lattice上求解沿着$x$或者$y$方向求解边界态的时候，还可以使用另外一种投影子的方法，这里就介绍一下这个方法。
{:.info}
<!--more-->
同样考虑BHZ模型加入Zeeman场这样的一个模型

$$
	H(\mathbf{k})=(m_0-t_x\cos k_x-t_y\cos k_y)\sigma_z+\lambda_x\sin k_x\sigma_xs_z+\lambda_y\sin k_y\sigma_y+h_0s_x
$$

这里$\Gamma_1=\sigma_zs_0,\Gamma_2=\sigma_xs_z,\Gamma_3=\sigma_ys_0,\Gamma_4=\sigma_0s_x$。
将$h_0s_{x}$视为微扰，求解BHZ模型对应的边界态。先沿着$y$方向取开放边界条件，进行Fourier变换到实空间

$$
H(k_{x})=(m_0-t_x\cos k_x-t_y/2(S_{y}+S_{y}^{\dagger}))\sigma_z+\lambda_x\sin k_x\sigma_xs_z+\lambda_y/(2i)(S_{y}-S_{y}^{\dagger})\sigma_y\label{eq8}
$$


$$
S_{y}\rvert i\rangle=\rvert i+1\rangle,\quad S_{y}^{\dagger}\rvert i\rangle=\rvert i-1\rangle
$$

考虑一个半无限大的系统，考虑$y$轴的正方向满足$S_y^\dagger\rvert 0\rangle=0$，此时可以见平移算符表示为

$$
S_y=\left[\begin{array}{ccccc}
0&0&0&0&\cdots\\
1&0&0&0&\cdots\\
0&1&0&0&\cdots\\
0&0&1&0&\cdots\\
\cdots&\cdots&\cdots&\cdots
\end{array}\right],\quad S^\dagger_y=\left[\begin{array}{ccccc}
0&1&0&0&\cdots\\
0&0&1&0&\cdots\\
0&0&0&1&\cdots\\
0&0&0&0&\cdots\\
\cdots&\cdots&\cdots&\cdots
\end{array}\right]
$$

考虑一个试探解$\rvert\psi_\mathbf{k}\rangle=\sum_{i=0}^{\infty}\lambda^i\rvert i\rangle\otimes\rvert\xi\rangle,\rvert\lambda\rvert<1$，求解(\ref{eq8})对应的本征方程，在$i\geq 1$的体态可以得到

$$
[(m_0-t_x\cos k_x-\frac{1}{2}(\lambda^{-1}+\lambda))\sigma_z+\lambda_x\sin k_x\sigma_xs_z+\frac{\lambda_y}{2i}(\lambda^{-1}-\lambda)\sigma_y]\rvert\xi\rangle=E\rvert\xi\rangle\label{eq9}
$$

在边界$i=0$的位置上，可以得到

$$
[(m_0-t_x\cos k_x-\frac{1}{2}\lambda^{-1})\sigma_z+\lambda_x\sin k_x\sigma_xs_z+\frac{\lambda_y}{2i}\lambda^{-1}\sigma_y]\rvert\xi\rangle=E\rvert\xi\rangle\label{eq10}
$$

通过联立(\ref{eq9})和(\ref{eq10})可以得到

$$
\sigma_x\rvert\xi\rangle=\rvert\xi\rangle\label{eq12}
$$

这表明边界态波函数的spinor部分对应的是$\sigma_xs_0$本征值为1的本征态，可以利用这些本征态构建投影算符

$$
\Pi^B=\frac{1}{2}(1+\sigma_xs_0)
$$

将投影算符作用在边界态的本征方程(\ref{eq10})上之后可以得到

$$
\lambda_x\sin k_x\sigma_xs_z\rvert\xi\rangle=E\rvert\xi\rangle\label{eq11}
$$

结合(\ref{eq10})和(\ref{eq11})以及(\ref{eq12})可以得到

$$
\lambda=m_0-t_x\cos k_x
$$

利用投影算符可以得到边界上的有效哈密顿量为

$$
\mathcal{H}_\text{eff}^B(\mathbf{k})=\Pi^BH(\mathbf{k})\Pi^B=\lambda_x\sin k_x\sigma_xs_z\Pi^B
$$

利用相同的方法，可以得到上边界的投影算符和有效哈密顿量

$$
\Pi^U=\frac{1}{2}(1-\sigma_zs_0),\quad\mathcal{H}_\text{eff}^U=\lambda_x\sin k_x\sigma_xs_z\Pi^U
$$


 当沿着$x$方向为开放边界的时候，可以分别求解得到左右两个边界上的投影算符和有效哈密顿量

$$
 \begin{aligned}
 \Pi^L&=\frac{1}{2}(1+\sigma_ys_z),\quad\mathcal{H}_\text{eff}^L=\lambda_y\sin k_y\sigma_y\Pi^L\\
 \Pi^R&=\frac{1}{2}(1-\sigma_ys_z),\quad\mathcal{H}_\text{eff}^R=\lambda_y\sin k_y\sigma_y\Pi^R
 \end{aligned}
$$

投影子作用在哈密顿量中的$\Gamma_i$矩阵上满足

$$
\begin{aligned}
&\Pi^{U/B}\Gamma_{1,3}\Pi^{U/B}=0,\quad\Pi^{R/L}\Gamma_{1,2}\Pi^{R/L}=0\\
&\Pi^{U/B}\Gamma_2\Pi^{U/B}=\Gamma_2\Pi^{U/B},\quad\Pi^{L/R}\Gamma_3\Pi^{L/R}=\Gamma_3\Pi^{L/R}\\
&\Pi^{U/B}\Gamma_4\Pi^{U/B}\neq0,\quad\Pi^{L/R}\Gamma_4\Pi^{L/R}=0
\end{aligned}
$$

根据上面的计算，当考虑$h_0\sigma_0s_x$层间耦合时，从上面的投影关系可以看到，这一项仅在$y$方向开边界的时候作为质量项存在，而对于$x$方向开边界的时候，这一项在投影算符的作用下为零。所以可以通过投影算子的方法得到两个方向上的边界有效哈密顿量为

$$
\begin{aligned}
H^{L/R}_{\text{eff}}&=\lambda_{y}\sin k_{y}\sigma_{y}\Pi^{L/R}\\
H^{U/B}_{\text{eff}}&=\lambda_{y}\sin k_{y}\sigma_{x}s_{z}\Pi^{U/B}+h_{0}\sigma_{x}\Pi^{L/R}\\
\end{aligned}
$$




即就是此时$k_x$的边界态被打开了一个能隙，但是$k_y$对应的边界态则始终是无能隙的，这与数值计算的结果是完全相同的。

# 参考
1.[Boundary Criticality of PT-Invariant Topology and Second-OrderNodal-Line Semimetals](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.126403)
