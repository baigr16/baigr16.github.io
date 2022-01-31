---
title: 哈密顿量的对称操作形式推导
tags: Study Group-Theory
layout: article
license: true
toc: true
key: a20211225
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
这里来详细的推导一下哈密顿量在对称操作下的形式以及理解一下为什么有时候要把对称操作表示为$C_mH(k)C_m^{-1}=H(R_mk)$的形式。
{:.info}
<!--more-->
# 哈密顿量对称变换
在看文章的过程中总是不太理解为什么要将哈密顿量满足的某种对称性表示为

$$
C_mH(k)C_m^{-1}=H(R_mk)
$$

的形式，这里就仔细推导并理解一下这个形式的由来。

首先系统的哈密顿量为$H$，它满足本征方程

$$
H\rvert\psi\rangle=E\rvert\psi\rangle
$$

假设现在有一个对称操作$P$要作用到系统上，那么就有

$$
PH\rvert\psi\rangle=EP\rvert\psi\rangle
$$

可以看到这个方程不再是一个本征方程了，而在量子力学中通常就喜欢处理本征值问题，所以将上面的形式改变一下

$$
PHP^{-1}\rvert\psi\rangle=E\rvert\psi\rangle
$$

那么这时候问题就变成了$PHP^{-1}$的本征方程，接下来的问题就是系统的哈密顿量在对称操作$P$下面到底是如何变化的。

在凝聚态系统里面，处理的哈密顿量都是带有动量$k$的$H(k)$，它与系统本来的哈密顿量H之间的关系是

$$
H(\mathbf{k})=e^{i\mathbf{k}\cdot\mathbf{R}}He^{-i\mathbf{k}\cdot\mathbf{R}}
$$

而$\hat{C}_m$是$H$的对称操作满足

$$
\hat{C}_mH\hat{C}_m^{-1}=H
$$

将这个操作作用到$H(\mathbf{k})$

$$
\hat{C}_me^{i\mathbf{k}\cdot\mathbf{R}}He^{-i\mathbf{k}\cdot\mathbf{R}}\hat{C}_m^{-1}=\hat{C}_me^{i\mathbf{k}\cdot\mathbf{R}}{\color{blue}\hat{C}_m^{-1}H\hat{C}_m}e^{-i\mathbf{k}\cdot\mathbf{R}}\hat{C}_m^{-1}
$$

下面的问题就是看看对称操作如何对$e^{-i\mathbf{k}\cdot\mathbf{R}}$作用，取一个本征态$\rvert\mathbf{r}\rangle$，有下面的关系

$$
\hat{C}_me^{i\mathbf{k}\cdot\mathbf{R}}\hat{C}_m^{-1}\rvert\mathbf{r}\rangle=\hat{C}_me^{i\mathbf{k}\cdot\mathbf{R}}\rvert R_m^{-1}\mathbf{r}\rangle=e^{i\mathbf{k}\cdot ({\color{blue}R_m^{-1}\cdot\mathbf{r}})}\hat{C}_m\rvert R_m^{-1}\mathbf{r}\rangle=e^{i\mathbf{k}\cdot ({\color{blue}R_m^{-1}\cdot\mathbf{r}})}\rvert\mathbf{r}\rangle
$$

$$
\hat{C}_me^{-i\mathbf{k}\cdot\mathbf{R}}\hat{C}_m^{-1}\rvert\mathbf{r}\rangle=\hat{C}_me^{-i\mathbf{k}\cdot\mathbf{R}}\rvert R_m^{-1}\mathbf{r}\rangle=e^{-i\mathbf{k}\cdot ({\color{blue}R_m^{-1}\cdot\mathbf{r}})}\hat{C}_m\rvert R_m^{-1}\mathbf{r}\rangle=e^{-i\mathbf{k}\cdot ({\color{blue}R_m^{-1}\cdot\mathbf{r}})}\rvert\mathbf{r}\rangle
$$

这里将矢量进行转动

$$
\mathbf{k}\cdot(R_m^{-1}\cdot\mathbf{r})=R_m\cdot\mathbf{k}\cdot R_m(R_m^{-1}\cdot\mathbf{r})=R_m\mathbf{k}\cdot\mathbf{r}
$$

这里的含义就是同时对两个矢量进行相同的转动操作，二者之间的点积是相同的。将上面这些关系综合起来就可以得到

$$
\hat{C}_me^{i\mathbf{k}\cdot\mathbf{R}}{\color{blue}\hat{C}_m^{-1}H\hat{C}_m}e^{-i\mathbf{k}\cdot\mathbf{R}}\hat{C}_m^{-1}=e^{iR_m\mathbf{k}\cdot\mathbf{R}}He^{-iR_m\mathbf{k}\cdot\mathbf{R}}=H(R_m\mathbf{k})=C_mH(\mathbf{k})C_m^{-1}
$$

所以这里最重要的就是推导出

$$
H(R_m\mathbf{k})=C_mH(\mathbf{k})C_m^{-1}\label{e1}
$$

这个关系表达式，这也算是解决了在看文章的时候一直很疑惑，为什么要将哈密顿量满足的一个对称操作表示成这样的形式。

# 对易关系
在得到了$H(R_m\mathbf{k})=C_mH(\mathbf{k})C_m^{-1}$这个关系时候，下面的问题就是它有什么用途，假如找到一个高对称的位置满足$R_m\mathbf{k}=\mathbf{k}$，那么将其带入到(\ref{e1})中可以有

$$
[\hat{C}_m,H(\mathbf{k})]=0
$$

在对称操作的高对称位置上，哈密顿量与对称操作算符是满足对易关系的，那么它们就会有共同本征态，从而就可以利用对称算符的性质来研究哈密顿量在高对称点的本征态的性质。

# 参考
- 1. 感谢崔先生和刘一辰的讨论。



