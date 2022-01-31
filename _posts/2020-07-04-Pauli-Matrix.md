---
title: Pauli矩阵及对称操作算符
tags: Study Topology
layout: article
license: true
toc: true
pageview: true
key: a202007041
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
在看文献的过程中通常都会遇到哈密顿量由Pauli矩阵来写出，而哈密顿量的对称操作同样可以通过Pauli矩阵来写出，这样可以直观的看出哈密顿量的对称性到底是什么样的，同时还可以根据对称性，利用Pauli矩阵之间的关系，搞清楚那些项是可以存在的，哪些项是不被对称性所允许的。
<!--more-->
# 基本关系
首先将[Pauli矩阵](https://en.wikipedia.org/wiki/Pauli_matrices)的基本关系罗列一下
$$\sigma_{1}^{2}=\sigma_{2}^{2}=\sigma_{3}^{2}=-i \sigma_{1} \sigma_{2} \sigma_{3}=\left(\begin{array}{cc}
1 & 0 \\
0 & 1
\end{array}\right)=I$$

$\det(\sigma_i)=0,tr(\sigma_i)=0$，每个Pauli矩阵的本征值都是$\pm1$，对应的本征矢量分别为

$$\begin{array}{ll}
\psi_{x+}=\frac{1}{\sqrt{2}}\left(\begin{array}{l}
1 \\
1
\end{array}\right), & \psi_{x-}=\frac{1}{\sqrt{2}}\left(\begin{array}{c}
1 \\
-1
\end{array}\right) \\
\psi_{y+1}=\frac{1}{\sqrt{2}}\left(\begin{array}{c}
1 \\
i
\end{array}\right), & \psi_{y-}=\frac{1}{\sqrt{2}}\left(\begin{array}{c}
1 \\
-i
\end{array}\right) \\
\psi_{z+}=\left(\begin{array}{c}
1 \\
0
\end{array}\right), & \psi_{z-}=\left(\begin{array}{c}
0 \\
1
\end{array}\right)
\end{array}$$

$\psi_{x+}=i\sigma_y\psi_{x-},\psi_{y+}=\sigma_x\psi_{y-},\psi_{z+}=\sigma_x\psi_{z-}$

# 对易及反对易关系

$[\sigma_a,\sigma_b]=2i\epsilon_{abc}\sigma_c,{\sigma_a,\sigma_b}=2\delta_{ab}I$

$$\begin{array}{ll}
{\left[\sigma_{1}, \sigma_{2}\right]=2 i \sigma_{3}} & \left\{\sigma_{1}, \sigma_{1}\right\}=2 I \\
{\left[\sigma_{2}, \sigma_{3}\right]=2 i \sigma_{1}} & \left\{\sigma_{1}, \sigma_{2}\right\}=0 \\
{\left[\sigma_{3}, \sigma_{1}\right]=2 i \sigma_{2}} & \left\{\sigma_{1}, \sigma_{3}\right\}=0 \\
{\left[\sigma_{1}, \sigma_{1}\right]=0} & \left\{\sigma_{2}, \sigma_{2}\right\}=2 I
\end{array}$$

**每个Pauli矩阵跟自己是对易的，跟其它两个分量的Pauli矩阵是反对易的**

这里再多说点，通常可能哈密顿量中有多个自由度，比如轨道($\sigma_i$)，自旋($s_i$)

当然了，不同自由度的Pauli矩阵之间肯定是没有关系的，不会存在上面所说的对易还是反对易的关系，因为本来它们就属于不同的自由度，所以自然也就不会去满足那些关系。代表同一自由度的Pauli矩阵之间才应该谈论对易和反对易的问题。

比如以这篇文章[Majorana Corner Modes in a High-Temperature Platform]( https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.096803 )为例，哈密顿量为$$H(\mathbf{k})= M(\mathbf{k}) \sigma_{z} \tau_{z}+A_{x} \sin k_{x} \sigma_{x} s_{z}+A_{y} \sin k_{y} \sigma_{y} \tau_{z} 
+\Delta(\mathbf{k}) s_{y} \tau_{y}-\mu \tau_{z}$$

它具有时间反演对成性$\mathcal{T} H(\mathbf{k}) \mathcal{T}^{-1}=H(-\mathbf{k})$，和粒子空穴对称性$\mathcal{C}H(\mathbf{k}) \mathcal{C}^{-1}=-H(-\mathbf{k})$，对称算符的具体形式为$\mathcal{T}=is_y\mathcal{K},\mathcal{C}=\tau_x\mathcal{K}$，这里的$\mathcal{K}$代表的是复共轭。下面来详细的验证一下这些对称性。

## [时间反演对称性]( https://en.wikipedia.org/wiki/T-symmetry )

简单的在哈密顿量上来说，时间反演要求动量相反$p\rightarrow-p,t\rightarrow-t$，自旋$\mathcal{T}S\mathcal{T}^{-1}=-S$，从算符的形式上看$\mathcal{T}c^\dagger_{ja\uparrow}\mathcal{T}^{-1}=c^\dagger_{ja\downarrow},\mathcal{T}c^\dagger_{ja\downarrow}\mathcal{T}^{-1}=-c^\dagger_{ja\uparrow},\mathcal{T}c^\dagger_{ka\sigma}\mathcal{T}^{-1}=i\sigma^y_{\sigma\sigma'}c_{-ka\sigma'}$(这里的$\sigma$暂时用来代表自旋的取向，不要和上面搞混)，这些就是时间反演基本性质的。当把哈密顿量用Pauli矩阵写出来后，基矢已经被抽离了出来。**有关哈密顿量的矩阵形式和算符形式关系可以参考[这里]( http://www.guanjihuan.com/archives/4867 )**。现在只需要验证这些对称性算符和哈密顿量之间的关系即可。

其实从时间反演的性质也可以很清楚的看明白为什么它的算符形式为$\mathcal{T}=is_y \mathcal{K}$

下面以时间反演对称性为例，验证一下哈密顿量的性质：

$$is_y\mathcal{K}*M({\mathbf{k}})\sigma_z\tau_z*(is_y\mathcal{K})^{-1}$$，从这里可以看到最后一个复共轭操作并没有具体的对象去作用，所以可以直接将它扔掉，而且它和矩阵求逆也没关系，而对$is_y$的求逆可以这样理解，对易Pauli矩阵$s_y$来说它是要求其逆矩阵，而对于虚数$i$来说，就是简单的求倒数而已，且第一个时间反演算符中的复共轭是要作用到它后面所的项的(**这里要强调一下，复共轭算符的作用是全局的，不在乎到底是哪个自由度，关心自由度的只有Pauli矩阵运算的时候，是同一自由度之间的Pauli矩阵进行运算**)，而$is_y$是个实数，所以它也相当于不进行任何操作(**因为$M({\mathbf{k}})\sigma_z\tau_z$同样也是实数**)，所以这个记过就简化为

$$M({\mathbf{k}})\sigma_z\tau_z*is_y*\frac{1}{i}s_y=M({\mathbf{k}})\sigma_z\tau_zs_ys_y=M({\mathbf{k}})\sigma_z\tau_z$$

而对于等式右边的结果是$M(\mathbf{-k})=M(\mathbf{k})$，所以这一项是符合时间反演不变的。$M(\mathbf{k})=m_0-2\sum_it_i\cos(k_i)$

$$is_y\mathcal{K}*A_x\sin(k_x)\sigma_xs_z*(is_y\mathcal{K})^{-1}=is_y\mathcal{K}*A_x\sin(k_x)\sigma_xs_z*\frac{1}{i}s_y=A_x\sin(k_x)\sigma_xis_ys_z\frac{1}{i}s_y$$

$$A_x\sin(k_x)\sigma_x*s_ys_zs_y=A_x\sin(k_x)\sigma_x*is_xs_y=A_x\sin(k_x)\sigma_x*i*is_z=-A_x\sin(k_x)\sigma_xs_z$$

这一项在等式的右边为$A_x\sin(-k_x)\sigma_xs_z=-A_x\sin(k_x)\sigma_xs_z$，结果和上面完全符合。

剩余的内容就不逐一验证了，类似的操作都是一样的。