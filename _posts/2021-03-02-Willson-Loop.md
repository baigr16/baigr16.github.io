---
title: Willson Loop学习笔记
tags: Topology
layout: article
license: true
toc: true
key: a20210302
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
这篇文章主要记录一下自己学习Willson loop时候的笔记,因为他实空间中的一些物理图像非常好用,而且也和拓扑不变量是有非常深刻的联系.
<!--more-->
虽然在前面我在计算拓扑不变量的时候,也使用过Willson loop,但是却没有仔细的对一些概念进行深入的理解,最近正好在看和高阶拓扑相关的文章,里面Willson loop起到了非常关键的作用,这里就整理一下自己的学习笔记.

在晶体中电子的分布是没有确定的位置的,都是以电子云来描述的,但是任然可以在确定的方向上电子出现的概率相对比较大,这里就先从实空间与栋梁空间两个角度来理解电子在一个周期的循环中,到底是如何变化的.

# 位置算符

首先实空间的位置算符可以表示为

$$\hat{x}=\sum_{R,\alpha}=c^\dagger_{R,\alpha}\rvert 0\rangle e^{-i\Delta_k(R+r_\alpha)}\langle 0\rvert c_{R,\alpha}$$

这里$\alpha$代表元胞中不同的轨道,$R$则代表着晶体中不同的元胞,$r_\alpha$表示元胞中不同轨道$\alpha$相对于正电荷中心的距离.$\Delta_k=2\pi/N$,$N$是晶体中元胞数目.
利用傅里叶变换

$$
\begin{equation}
\begin{aligned}
c_{R,\alpha}&=\frac{1}{\sqrt{N}}\sum_ke^{-ik(R+r_\alpha)}c_{k,\alpha},\\
c_{k,\alpha}&=\frac{1}{\sqrt{N}}\sum_Re^{ik(R+r_\alpha)}c_{r,\alpha}
\end{aligned}
\end{equation}
$$

这里的$k\in \Delta_k\cdot(0,1,\cdots,N-1)$.采用周期性边界条件

$$c_{R+N,\alpha}=c_{R,\alpha}\rightarrow c_{k+G,\alpha}=e^{iGr_\alpha}c_{k,\alpha}$$

这里的$G$是倒空间的晶格矢量.在这个新的基矢下可以将位置算符改写为

$$\hat{x}=\sum_{k,\alpha}c^\dagger_{k+\Delta_k,\alpha}\rvert 0\rangle\langle 0\rvert c_{k,\alpha}$$

同样的,二次量子化的哈密顿量可以写成

$$H=\sum_k c^\dagger_{k,\alpha}\left[h_k\right]^{\alpha,\beta}c_{k,\beta}\label{p2}$$

因为前面采用了周期边界条件,哈密顿量满足

$$h_{k+G}=V^{-1}(G)h_kV(G)\qquad \left[V(G)\right]^{\alpha,\beta}=e^{-iGr_\alpha}\delta_{\alpha,beta}\label{p1}$$

将哈密顿量进行对角化

$$\left[h_k\right]^{\alpha,\beta}=\sum_n\left[u^n_k\right]^\alpha\epsilon_{n,k}\left[u^{*n}_k\right]^\beta$$

这里的$\left[u^n_k\right]^\alpha$是本征矢量$\rvert u_k^n\rangle$的第$\alpha$个分量.

为了强制满足周期边界条件(\ref{p1}),这里采用周期规范

$$\left[u_{k+G}^n\right]^\alpha=\left[V^{-1}(G)\right]^{\alpha,\beta}\left[u^n_k\right]^\beta$$

这时候可以将哈密顿量(\ref{p2})表示为

$$\begin{equation}
\begin{aligned}
H&=\sum_{\alpha,k}\gamma^\dagger_{n,k}\epsilon_{n,k}\gamma_{n,k}\\
\gamma_{n,k}&=\sum_\alpha\left[u^{*n}_k\right]^\alpha c_{k,\alpha}
\end{aligned}
\end{equation}\label{p3}$$

由于是在周期条件下,算符$\gamma$满足

$$\gamma_{n,k}=\gamma_{n,k+G}$$

**因为我们这里关心的主要是绝缘体在零温下砸情况,所以只关心占据态的电子能带.**构建占据态能带的投影算符

$$P^{occ}=\sum_{n=1}^{N_{occ}}\sum_k\gamma_{n,k}^{\dagger}\rvert 0\rangle\langle 0\rvert\gamma_{n,k}$$

**接下来就将位置算符投影到占据态子空间,然后进行对角化**

$$P^{occ}\hat{x}P^{occ}=\sum_{n,k}\sum_{n^{'},k^{'}}\gamma^\dagger_{n,k}\rvert 0\rangle\langle 0\rvert\gamma_{n^{'}k^{'}}(\sum_{q,\alpha}\langle 0\rvert\gamma_{m,k}c^\dagger_{q+\Delta_k,\alpha}\rvert 0\rangle\langle 0\rvert c_{q,\alpha}\gamma^\dagger_{n^{'}k^{'}}\rvert 0\rangle)$$

从公式(\ref{p3})可以得到$\langle 0\rvert\gamma_{n,k}c^\dagger_{q,\alpha}\rvert 0\rangle=\left[u^{*n}_k\right]^\alpha\delta_{k,q}$,可以将投影后的位置算符约化成

$$P^{occ}\hat{x}P^{occ}=\sum_{m,n=1}^{N_{occ}}\sum_k\gamma^\dagger_{m,k+\Delta_k}\rvert 0\rangle\langle u^m_{k+\Delta_k}\rvert u^n_k\rangle\langle 0\rvert\gamma_{n,k}$$

这里的一些符号简写为$\langle u^m_q\rvert u^n_k\rangle=\sum_\alpha\left[u^{*m}_q\right]^\alpha\left[u^n_k\right]^\alpha$

矩阵$G$的分量形式为$\left[G_k\right]^{mn}=\langle u^m_{k+\Delta_k}\rvert u^n_k\rangle$,接下来对矩阵$G$进行奇异值分解

$$G=UDV^\dagger$$

这里的矩阵$D$只有对角线元素不为0.**The failure of G to be unitary
is manifest in the fact that the (real-valued) singular values
along the diagonal of D are less than 1  **,所以对于每一个$k$可以定义

$$F=UV^\dagger$$

**在热力学极限在$N\rightarrow\infty$**时有$\left[F\right]^{mn}=\left[G\right]^{mn}$.

为了对角化投影位置算符,可以写出下列本征值问题

$$(P^{occ}\hat{x}P^{occ})\rvert\Psi^j\rangle=E^j\rvert\Psi^j\rangle$$

在$\gamma_{n,k}\rvert 0\rangle$这个基矢下,其对应的矩阵形式为

$$\left(\begin{array}{ccccc}0&0&0&\cdots&F_{k_N}\\F_{k_1}&0&0&\cdots&0\\0&F_{k_2}&0&\cdots&0\\\cdots&\cdots&\cdots&\cdots&\cdots\\0&0&0&\cdots&0\end{array}\right)\left(\begin{array}{c}\nu_{k_1}\\\nu_{k_2}\\\nu_{k_3}\\\cdots\\\nu_{k_N}\end{array}\right)^j=E^j\left(\begin{array}{c}\nu_{k_1}\\\nu_{k_2}\\\nu_{k_3}\\\cdots\\\nu_{k_N}\end{array}\right)^j$$

这里$k_N=\Delta_k(N-1),j\in 1\cdots N_{occ}$.在这里利用$F_k$代替了$G_k$是为了保证幺正性.通过重复的利用上面的方程,我们可以得到下面的关系

$$\mathcal{W}_{k_f\leftarrow k_i}\rvert\nu^j_{k_i}\rangle=(E^j)^{(k_f-k_i)/\Delta_k}\rvert\nu^j_{k_f}\rangle\label{p4}$$

定义离散的Willson line

$$\mathcal{W}_{k_f\leftarrow k_i}=F_{k_f-\Delta_k}F_{k_f-\Delta_k}F_{k_f-2\Delta_k}\cdots F_{k_i+\Delta_k}F_{k_i}$$

可以让这个Willson line穿过整个BZ,则(\ref{p4})对应的本征值问题就变为

$$\mathcal{W}_{k+2\pi\leftarrow k}\rvert v_k^j\rangle=(E^j)^N\rvert v_k^j\rangle$$

**由于Willson loop是幺正的,所以它的本征值是个简单的位相**

$$(E^j)^N=e^{i2\pi v^j}$$

它有$N$个解

$$E^{j,R}=e^{i2\pi v^j/N+i2\pi R/N}=e^{i\Delta_k(v^j+R)}$$

**这里的位相$v^j$就是Wannier Center,代表着在一个元胞中电子相对于电荷中心偏离的位置**

# 极化与拓扑不变量的关系

The electronic contribution to the dipole moment,measured as the electron charge times the displacement of the electrons from
the center of the unit cell

$$p=\sum_j v^j$$

利用Willson loop可以将电子的极化写成

$$p=-\frac{i}{2\pi}\textrm{log det}\left[\mathcal{W}_{k+2\pi\leftarrow k}\right]$$

Berry曲率

$$\left[\mathcal{A}_k\right]^{mn}=-i\langle u_k^m\rvert\partial_k\rvert u_k^n\rangle$$

$$p=-\frac{i}{2\pi}\textrm{log det}\left[e^{-i\int_k^{k+2\pi}\mathcal{A}_kdk}\right]=-\frac{1}{2\pi}\int_k^{k+2\pi}\textrm{Tr}\left[\mathcal{A}_k\right]dk\qquad \textrm{mod}\quad 1$$

这就是Willson loop与拓扑不变量之间的联系.

# 参考文献
1.[Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)


# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)