---
title: Topological insulators with inversion symmetry
tags: Topology
layout: article
license: true
toc: true
key: a20210506
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
这里整理一下当拓扑绝缘体存在时间反演对称之后,$\mathcal{Z}_2$拓扑不变量可以只通过在时间反演不变动量点$\Gamma_i$计算宇称,这极大的简化了$\mathcal{Z}_2$的计算,而且在材料计算中也有很重要的作用.
<!--more-->

首先通过[时间反演极化](https://yxli8023.github.io/2021/04/25/TR-Polarization.html)这篇博客,先了解了时间反演不变的系统中,可以通过定义$\mathcal{Z}_2$拓扑不变量来表征系统的拓扑性质,但是真实的计算$\mathcal{Z}_2$的过程中会涉及到复杂的规范选择问题,而且通常这个规范是不容易找到的,如果系统存在**反演对称**时,$\mathcal{Z}_2$的计算就会得到极大的简化.

# $\mathcal{Z}_2$ invariant with inversion symmetry
假设系统存在一个反演中心${\bf r}=0$,且满足反演对称性

$$[\mathcal{H},P]=0\qquad H(-{\bf k})=PH({\bf k})P^{-1}\qquad P\rvert {\bf r},s_z\rangle=\rvert {\bf -r},s_z\rangle$$

这里${\bf r}$是三维空间坐标,$s_z$表示自旋,**因为自旋是赝矢量,所以在空间反演下是不会改变的**.

当系统存在时间反演时贝利矢势满足$\mathcal{F}({\bf -k})=-\mathcal{F}({\bf k})$,存在空间反演时$\mathcal{F}({\bf -k})=\mathcal{F}({\bf k})$,贝利曲率为$\nabla_k\times \mathcal{A}({\bf k})$

$$\mathcal{A}({\bf k})=-i\sum_{n=1}^{2N}\langle u_{n,{\bf k}}\rvert\nabla_k\rvert u_{n,{\bf k}}\rangle$$

这里需要对所有的占据态求和$2N$,从上面可以知道,当同时存在空间反演与时间反演时$\mathcal{F}({\bf k})=0$,则可以选择一个全局连续的"横场"规范$\mathcal{A}({\bf k})=0$.

在任意规范下,考虑一个$2N\times 2N$的矩阵

$$v_{mn}({\bf k})=\langle u_{m,{\bf k}}\rvert P\Theta\rvert u_{n,{\bf k}}\rangle$$

由于$\langle a\rvert b\rangle=\langle\Theta b\rvert\Theta a\rangle$以及$\Theta^2=-1$,所以矩阵$v({\bf k})$是个反对称矩阵,而且$[P\Theta,H({\bf k})]=0$,所以$v({\bf k})$是幺正的,因此可以定义$v({\bf k})$的Pfaffian,并且它是幺模的.$\text{Pf}[v({\bf k})]$的位相是和规范选择相关的,它的梯度和$\mathcal{A}({\bf k})$是紧密相关的.

$$\mathcal{A}({\bf k})=-\frac{i}{2}\text{Tr}[v({\bf k})^\dagger\nabla_kv({\bf k})]=-i\nabla_k\log[\text{Pf}[v({\bf k})]]$$

上面的推导中用到了$\qquad\text{det}[v]=\text{Pf}[v]^2\qquad\nabla_k\log[\text{det}[v]]=\text{Tr}[\nabla_k\log[v({\bf k})]]=\text{Tr}[v^\dagger({\bf k})\nabla_kv({\bf k})]$

为了满足$\mathcal{A}({\bf k})=0$这个条件,只需要调整$\rvert u_{n\bf k}\rangle$的位相,使其满足

$$\text{Pf}[v({\bf k})]=1$$

在计算$\mathcal{Z}_2$的时候,需要计算

$$w_{mn}({\bf k})\equiv\langle u_{m\bf -k}\rvert\Theta\rvert u_{n\bf k}\rangle$$

再结合空间反演操作$P$,可得到

$$w_{mn}(\Gamma_i)=\langle\psi_{m,\Gamma_i}\rvert P(P\Theta)\rvert\psi_{n,\Gamma_i}\rangle$$

这里$\Gamma_i$表示时间反演不变动量点,$P^2=1$,将$\rvert u_{n\Gamma_i}\rangle$代替为$\rvert\psi_{n\Gamma_i}\rangle=\rvert\psi_{n-\Gamma_i}\rangle$.因为$[\mathcal{H},P]=0,\rvert\psi_{n\Gamma_i}\rangle$是$P$的本征态,对应的本征值为$\xi_n(\Gamma_i)=\pm 1$,当把$\rvert \psi_{n\Gamma_i}\rangle$回代为$\rvert u_{n\Gamma_i}\rangle$之后

$$w_{mn}(\Gamma_i)=\xi_m(\Gamma_i)v_{mn}(\Gamma_i)$$

矩阵$w$的Pfaffian满足

$$\text{Pf}[w]^2=\text{det}[w]=\text{det}[v]\Pi_{n=1}^{2N}\xi_n\label{eq1}$$

由于Kiamers简并的存在$\rvert u_{2m,\Gamma_i}\rangle$与$\rvert u_{2m+1,\Gamma_i}\rangle=\Theta\rvert u_{2m,\Gamma_i}\rangle$具有相同的宇称本征值,则(\ref{eq1})的乘积中,每个本征值都会出现两次,取平方根可得

$$\text{Pf}[w]=\text{Pf}[v]\Pi_{m=1}^{N}\xi_{2m}$$

由于$\text{Pf}[v]=1$,则在"横场"规范下,可以得到

$$\delta_i=\frac{\sqrt{\text{det}[w(\Gamma_i)]}}{\text{Pf}[w(\Gamma_i)]}=\pm 1=\Pi_{m=1}^{N}\xi_{2m}(\Gamma_i)$$

对于二维系统,单个$\mathcal{Z}_2$不变量表示为

$$(-1)^\nu=\Pi_{i=1}^4\delta_i$$

在三维有8个时间反演不变动量点,将会存在4个独立的$\mathcal{Z}_2$不变量,其中的$\nu_0$表示为8个时间反演不变动量点的乘积

$$(-1)^\nu=\Pi_{i=1}^8\delta_i$$

其余的三个则是处在相同平面上4个时间反演不变动量点的乘积

$$(-1)^{\nu_k}=\Pi_{n_k=1;n_{j\neq k=0,1}}\delta_{i=(n_1,n_2,n_3)}$$

# 参考
- 1.[Topological insulators with inversion symmetry](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.045302)