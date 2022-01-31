---
title: cylinder结构上的边界态理论
tags: Topology 
layout: article
license: true
toc: true
key: a20211215
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
这里整理一下如何在一个方向上开边界的时候，如何求解边界态理论。
{:.info}
<!--more-->
# 前言
之前已经整理过在动量空间中求解边界态理论，这里整理一下利用另外的一种方法来求解边界态理论，主要就是将哈密顿量变换到cylinder结构上，在一个有效大小的方向上计算边界态理论，主要的参考[有效边界理论(space部分)](https://yxli8023.github.io/2021/01/20/Effective-Edge-Theory.html)这篇Blog中的哈密顿量来进行研究。文章中的主要内容都是参考[Majorana Corner Modes in a High-Temperature Platform](https://arxiv.org/abs/1803.08545)这篇文章中的内容。
# 模型
这里采用文章中模型来进行求解

$$H(\mathbf{k})=(m_0-t_x\cos k_x-t_y\cos k_y)\sigma_z\tau_z+A_x\sin k_x\sigma_xs_z+A_y\sin k_y\sigma_y\tau_z+[\Delta_0-\Delta_1(\cos k_x+\cos k_y)]s_y\tau_y$$

首先沿着一个方向进行傅立叶变换，使得一个方向在动量空间另外一个方向在实空间，将哈密顿量表示为

$$
\begin{equation}
\begin{aligned}
H=&\sum_{x,k_y}\Psi^\dagger_{x,k_y}[(m_0-t_y\cos k_y)\sigma_z\tau_z+A_y\sin k_y\sigma_y\tau_z+(\Delta_0-\Delta_1\cos k_y)s_y\tau_y]\Psi_{x,k_y}\\
&+\sum_{x,k_y}[\Psi^\dagger_{x,k_y}(-\frac{t_x}{2}\sigma_z\tau_z+i\frac{A_x}{2}\sigma_xs_z-\frac{\Delta_1}{2}s_y\tau_y)\Psi_{x+1,k_y}+\text{h.c}]
\end{aligned}
\end{equation}
$$

这里基矢的选择为

$$
\Psi_{x,k_y}=(c_{x,a,k_y,\uparrow},c_{x,b,k_y,\uparrow},c_{x,a,k_y,\downarrow},c_{x,b,k_y,\downarrow},c^\dagger_{x,a,-k_y,\uparrow},c^\dagger_{x,b,-k_y,\uparrow},c^\dagger_{x,a,-k_y,\downarrow},c^\dagger_{x,b,-k_y,\downarrow})
$$

这的$x$表示的就是沿着$x$方向的格点的索引，这里将哈密顿量的基矢选择为$(\Psi_{1,k_y},\Psi_{2,k_y},\cdots)$，可以将哈密顿量表示为$H(k_y)=H_0(k_y)+H_1(k_y)+H_2(k_y)$，可以将其表示为


$$
\begin{equation}
H_0(k_y)=\left[
\begin{array}{cccc}
D_0&T_0&0&\cdots\\
T^\dagger_0&D_0&T_0&\cdots\\
0&T^\dagger_0&D_0&\cdots\\
\cdots&\cdots&\cdots&\cdots
\end{array}
\right]
\end{equation}\quad D_0=(m_0-t_y\cos k_y)\sigma_z\tau_z,T_0=(-t_x\sigma_z\tau_z+iA_x\sigma_xs_z)/2
$$

$$
\begin{equation}
H_1(k_y)=\left[
\begin{array}{cccc}
D_1&T_1&0&\cdots\\
T^\dagger_1&D_1&T_1&\cdots\\
0&T^\dagger_1&D_1&\cdots\\
\cdots&\cdots&\cdots&\cdots
\end{array}
\right]
\end{equation}\quad D_1=(\Delta_0-\Delta_1\cos k_y)s_y\tau_y,T_1=-\Delta_1s_y\tau_y/2
$$

$$
\begin{equation}
H_2(k_y)=\left[
\begin{array}{cccc}
A_y\sin k_y\sigma_y\tau_z&0&0&\cdots\\
0&A_y\sin k_y\sigma_y\tau_z&0&\cdots\\
0&0&A_y\sin k_y\sigma_y\tau_z&\cdots\\
\cdots&\cdots&\cdots&\cdots
\end{array}
\right]
\end{equation}
$$
同样的，在这里将$H_1(k_y),H_2(k_y)$看作时微扰，而首先来求解$H_0(k_y=0)$的本征态，可以将其本征波函数表示为

$$
\psi=(\psi_1,\psi_2,\cdots)
$$

将哈密顿量代入即可以得到其满足迭代方程

$$
\frac{1}{2}(-t\sigma_z\tau_z+iA\sigma_xs_z)\psi_{n+1}+m\sigma_z\tau_z\psi_n+\frac{1}{2}(-t\sigma_z\tau_z-iA\sigma_xs_z)\psi_{n-1}=0\label{f1}
$$

这里做了参数简化$m_0-t_x\rightarrow m,t_x\rightarrow t,A_x\rightarrow A$，同样的可以将空间部分和spinor部分独立求解，这里可以看到有下面的关系

$$
{\sigma_z\tau_z,\sigma_xs_z}=0\quad \Gamma_1=\sigma_z\tau_z,\Gamma_2=\sigma_xs_z
$$

从而可以选择一个平庸的解

$$
\psi_{n,\alpha}=\lambda^n\phi_\alpha
$$

这里波函数的spinor部分满足

$$
-i\Gamma_1\Gamma_2\phi_\alpha=s\phi_\alpha\quad s=\pm 1
$$

将这个波函数代入(\ref{f1})中就可以得到

$$
\frac{1}{2}(-t+As)\lambda^{n+1}+m\lambda^n+\frac{1}{2}(-t-As)\lambda^{n-1}=0
$$

为了让方程具有非平庸解，需要令$\lambda\neq0$，从而可以将方程化简为一个二次方程


$$
\frac{1}{2}(-t+As)\lambda^2+m\lambda+\frac{1}{2}(-t-As)=0
$$

可以求解得到$\lambda$的表达式为

$$
\lambda_\pm=\frac{m\pm\sqrt{m^2-(t^2-A^2)}}{t+As}
$$

当考虑一个半无限大的系统的时候$L\rightarrow\infty$，在边界上其波函数满足

$$
\psi_{\alpha,0}=\psi_{\alpha,+\infty}=0
$$

这里人为的在$x=0$的位置上加入了一个格点位置，所以可以将解的形式表示为

$$
\psi_{\alpha,n}=\mathcal{N}(\lambda^n_+-\lambda_{-}^n)\phi_\alpha
$$

因为波函数要满足归一化，所以

$$
\psi^\dagger*\psi=1\quad \psi=(\psi_1,\psi_2,\cdots)^T
$$

可以得到

$$
\mathcal{N}^2\sum_{n=0}^\infty(\lambda_{+}-\lambda_{-})^2=1
$$

从而可以得到归一化系数为

$$
\rvert\mathcal{N}\rvert^2=[\frac{ts}{A}\frac{\rvert m^2-(t^2-A^2)}{t^2-m^2}]^{-1}
$$

得到了这些结果之后，就可以进行微扰计算了，将$H_1,H_2$投影到$H_0$的本征态上，并且忽略$k_y$以上的项。

$$
H_{1,\alpha\beta}=\psi^\dagger_\alpha[H_1(k_y)+H_2(k_y)]\psi_\beta
$$

可以对微扰的哈密顿量进行分析可以知道，其实也只有$H_1(k_y)$项才会产生质量项，而且在这里需要忽略除线性项之外的其它项，所以在投影的时候得到的质量项为

$$
M_1=\Delta_0-\Delta_1-\frac{\Delta_1}{2}(\sum_{n=1}\psi^\dagger_{n+1}\psi_n+\sum_{n=2}\psi^\dagger_{n-1}\psi_n)
$$

这里的

$$
(\sum_{n=1}\psi^\dagger_{n+1}\psi_n+\sum_{n=2}\psi^\dagger_{n-1}\psi_n)=\mathcal{N}^2[(\sum_{n=1}(\lambda^{n+1}_+-\lambda^{n+1}_{-})(\lambda^{n}_+-\lambda^{n}_{-})+\sum_{n=2}(\lambda^{n-1}_+-\lambda^{n-1}_{-})(\lambda^{n}_+-\lambda^{n}_{-})]=\mathcal{N}^2\frac{ms}{A}\frac{\rvert m^2-(t^2-m^2) \rvert}{t^2-m^2}
$$

代入参数就可以得到

$$
\mathcal{N}^2\frac{ms}{A}\frac{\rvert m^2-(t^2-m^2) \rvert}{t^2-m^2}=\frac{m}{t}
$$

将原本的参数简化代入就可以得到

$$
M_1=\Delta_0-\Delta_1-\Delta_1\frac{m_0-t_x}{t_x}
$$

这个结果与之前通过动量空间中的哈密顿量通过边界态理论得到的结果是相同的。这里在级数求和的时候用了一下软件计算，结算代码如下图所示

![png](/assets/images/Mma/lattice-edge-theory.png)

# 参考
- 1.[Majorana Corner Modes in a High-Temperature Platform](https://arxiv.org/abs/1803.08545)

