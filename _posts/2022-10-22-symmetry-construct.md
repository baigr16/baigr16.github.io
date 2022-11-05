---
title: 根据哈密顿量获得其对称操作算符
tags:  Mathematica Topolog
layout: article
license: true
toc: true
key: a20221022
pageview: true
# cover: /assets/images/GroupTheory/cube_symmetry.jpg
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#203028'
  background_image:
    gradient: 'linear-gradient(135deg, rgba(34, 139, 87 , .4), rgba(139, 34, 139, .4))'
aside:
  toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
通常给出的一个哈密顿量，它的对称性可能是比较高的，有时候在看文章的时候可能也并不会将哈密顿量满足的所有对称性都列举出来，但是在分析问题的时候总是会用到这些对称性，这里就整理个程序分析一下如果哈密顿量具有某种对称性，那么它对应的矩阵形式是怎样的。其实这里给出的也不一定就是正确的，因为我并不知道哈密顿量的基矢是什么，能给出的也就是一些操作矩阵满足对称操作对哈密顿量的变换。
{:.info}
<!--more-->
# 简略分析
假如我现在有一个哈密顿量

$$\begin{equation}H(\mathbf{k})=(m_0-t_x\cos k_x-t_y\cos k_y)\sigma_z+A_x\sin k_x\sigma_xs_z+A_y\sin k_y\sigma_y +w_0(\cos k_x-\cos k_y)\sigma_xs_x\label{q1}\end{equation}$$

选择一组参数

$$m_0=1.5,t_x=t_y=A_x=A_y=1.0, w_0=0.5$$

这里我也先不明确哈密顿量中的不同Pauli算符代表的自由度是什么，因为我只想找到满足对称操作的矩阵，可能得到的并不是真正的操作算符。假如我想找哈密顿量\eqref{q1}的镜面操作$\mathcal{M}_x$和$\mathcal{M}_y$，我先不关心它是不是真的存在这两个对称性，反正它们对哈密顿量的操作满足

$$\begin{align}&\mathcal{M}_xH(k_x,k_y)\mathcal{M}_x^{-1}=H(-k_x,k_y)\\ &\mathcal{M}_yH(k_x,k_y)\mathcal{M}_y^{-1}=H(k_x,-k_y)\\\end{align}$$

反正这个操作算符的维度和哈密顿量的维度肯定是相同的，那么就将所有16个Pauli矩阵的直积全部构建出来。但是实际上应该是有32个，因为如果我们考虑的是spinfull的系统，那么此时$\mathcal{M}_i^2=-1$，此时还需要考虑16个Pauli矩阵的直积再乘以一个虚数$i$，那么就先将这32个Pauli矩阵作为操作元，分别作用到哈密顿量的每一项上面。

比如

$$(m_0-t_x\cos k_x-t_y\cos k_y)\rightarrow (m_0-t_x\cos (-k_x)-t_y\cos(k_y))\rightarrow (m_0-t_x\cos k_x-t_y\cos k_y)$$

在镜面$\mathcal{M}_x$下面系数不变号，那么矩阵$\mathcal{M}_x$对系数后面的Pauli矩阵的操作即满足

$$\mathcal{M}_x(\sigma_zs_0)\mathcal{M}_x^{-1}=\sigma_zs_0$$

其余项的分析也是完全相同的

$$\begin{align}&\sin(k_x)\rightarrow\sin(-k_x)\rightarrow -\sin(k_x)\quad \mathcal{M}_x(\sigma_xs_z)\mathcal{M}_x^{-1}=-\sigma_xs_z \\ &\sin(k_y)\rightarrow\sin(k_y)\rightarrow \sin(k_y)\quad \mathcal{M}_x(\sigma_ys_0)\mathcal{M}_x^{-1}=\sigma_ys_0\\ &\cos(k_x)\rightarrow \cos(-k_x)\rightarrow\cos(k_x)\quad \mathcal{M}_x(\sigma_xs_x)\mathcal{M}_x^{-1}=\sigma_xs_x\end{align}$$

通过上面的这些对易以及反对易关系，就可以确定$\mathcal{M}_x$的具体形式，但是需要注意的是，因为我们可能处理的是spinless或者spinfull的体系，所以到底镜面对称操作中是否存在虚数$i$就要根据具体情况来选择了。我这里给出的程序同时包含了这两种情况，所以对于具体问题的分析，还是需要具体对待，这里也就是给了一点分析的方法。

# 代码
因为无法在Blog中粘贴Mathematica的代码，所以截图示意一下，完整的代码<a class="button button--success button--rounded button--lg" href="/assets/data/symmetry.nb"><i class="fas fa-download"></i> 点击这里下载</a>


![png](/assets/images/Mma/symmetry.png)

# 补充
在[根据对称性计算体系电四极矩](https://yxli8023.github.io/2022/10/21/symmetry-quadrupole.html)这篇Blog中，就是通过这个方法，对BHZ+$d$波配对的模型构建得到了镜面对称操作算符，而且发现[Detection of second-order topological superconductors by Josephson junctions](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.012018)这篇文章中给出的镜面操作算符其实是错的，因为这个模型应该是spinfull的，所以镜面操作满足

$$\mathcal{M}_i^2=-1$$

而且通过计算电四极矩也发现，它文章中给出的镜面对称算符根本得到不$Q_{xy}$的量子化值，得到的始终是零。