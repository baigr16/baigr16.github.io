---
title: 重学Wilson loop
tags: Topology 
layout: article
license: true
toc: true
key: a20221010
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
这里整理一下自己在重新学习高阶拓扑过程中，对Wilson loop以及相关内容的一个重新理解，也为后面理解Nested Wilson loop做一个铺垫。
{:.info}
<!--more-->

# 重学Wilson loop
最近重读经典文献[Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)，对Wilson loop其中的一些表示还有具体含义有了更多的理解，之前都只是参照公式进行计算，现在想重新整理一下具体的物理含义，而不是简简单单的计算一下。
# 位置算符
在理解Berry位相的时候，除了通过参数空间中的演化来理解，在晶体中的电极化也是跟Berry位相有关的，实空间中极化问题的解决也就是现代极化理论。其中关键的一定就是将位置算符在周期系统中表示出来，首先在实空间中写出了位置算符，假设晶体中有$N$个原胞，每个原胞中的轨道数量为$N_{\rm orb}$，可以将位置算符表示为

$$\hat{x}=\sum_{R,\alpha}c^\dagger_{r,\alpha}\rvert 0\rangle e^{-i\delta_k(R+r_\alpha)}\langle 0\rvert c_{R,\alpha}$$

这里的$\alpha\in1\cdots N_{\rm orb}$是轨道的指标，$R$则是原胞索引，$r_\alpha$是原胞中原子相对于原胞中心的位置，$\Delta_k=2\pi/N$，想要研究其在动量空间中的性质，进行傅立叶变换即可

$$\begin{array}{l}c_{R, \alpha}=\frac{1}{\sqrt{N}} \sum_{k} e^{-i k\left(R+r_{\alpha}\right)} c_{k, \alpha}, \\c_{k, \alpha}=\frac{1}{\sqrt{N}} \sum_{R} e^{i k\left(R+r_{\alpha}\right)} c_{R, \alpha}\end{array}$$

这里的$k\in\Delta_k\cdot (0,1,\cdots,N-1)$，变换过程中利用周期边界条件

$$c_{R+N, \alpha}=c_{R, \alpha} \rightarrow c_{k+G, \alpha}=e^{i G r_{\alpha}} c_{k, \alpha}$$

在动量空间中就可以将位置算符表示为

$$\hat{x}=\sum_{k, \alpha} c_{k+\Delta_{k}, \alpha}^{\dagger}|0\rangle\langle 0| c_{k, \alpha}$$

在二次量子化形式下，哈密顿量表示为

$$H=\sum_{k} c_{k, \alpha}^{\dagger}\left[h_{k}\right]^{\alpha, \beta} c_{k, \beta},$$

通常情况下我们在计算的时候，输入其实是中间的矩阵部分$\left[h_{k}\right]^{\alpha, \beta}$，而且在求解本征值以及本征态的时候，也是在对角化这个矩阵，可以将对角化之后的哈密顿量表示为

$$\left[h_{k}\right]^{\alpha, \beta}=\sum_{n}\left[u_{k}^{n}\right]^{\alpha} \epsilon_{n, k}\left[u_{k}^{* n}\right]^{\beta},$$

此时二次量子化形式的哈密顿量为

$$H=\sum_{n, k} \gamma_{n, k}^{\dagger} \epsilon_{n, k} \gamma_{n, k}$$

这里的$\gamma_{n,k}$就是对角表象下面的准粒子，它与电子产生算符之间的变换关系为

$$\gamma_{n, k}=\sum_{\alpha}\left[u_{k}^{* n}\right]^{\alpha} c_{k, \alpha}$$

其实就是通过本征矢量将电子算符表象和准粒子$\gamma_{n,k}$表象联系起来。

# 占据态空间
在研究拓扑问题的时候，通常要关心的都是占据态的性质。在前面通过对角化的方式得到了系统所有的本征态，那么就可以利用所有的占据态构建投影算符

$$P^{\text {occ }}=\sum_{n=1}^{N_{\text {occ }}} \sum_{k} \gamma_{n, k}^{\dagger}|0\rangle\langle 0| \gamma_{n, k},$$

这里$N_{\rm occ}$就是占据态的数量。接下来就是要将位置算符投影到占据态空间中

$$P^{\mathrm{occ}} \hat{x}P^{\mathrm{occ}}= \sum_{n, k} \sum_{n^{\prime}, k^{\prime}} \gamma_{n, k}^{\dagger}|0\rangle\langle 0| \gamma_{n^{\prime}, k^{\prime}}\left(\sum_{q, \alpha}\left\langle 0\left|\gamma_{n, k} c_{q+\Delta_{k}, \alpha}^{\dagger}\right| 0\right\rangle\left\langle 0\left|c_{q, \alpha} \gamma_{n^{\prime}, k^{\prime}}^{\dagger}\right| 0\right\rangle\right)$$

现在就需要利用准粒子算符和电子算符之间的变换关系了$\gamma_{n, k}=\sum_{\alpha}\left[u_{k}^{* n}\right]^{\alpha} c_{k, \alpha}$，利用这个关系可以得到

$$\langle 0\rvert \gamma_{n,k}c^\dagger_{q,\alpha}\rvert 0\rangle=[u_k^{*n}]^{\alpha}\delta_{k,q}$$

利用这个关系可以将位置算符化简为

$$P^{\mathrm{occ}} \hat{x} P^{\mathrm{occ}}=\sum_{m, n=1}^{N_{\text {occ }}} \sum_{k} \gamma_{m, k+\Delta_{k}}^{\dagger}\rvert 0\rangle\left\langle u_{k+\Delta_{k}}^{m} \mid u_{k}^{n}\right\rangle\langle 0\rvert \gamma_{n, k},$$

在上面的表达式中，会有一点让人误解的地方，其实

$$\langle u^m_q\rvert u^n_k\rangle=\sum_\alpha[u_q^{*m}]^\alpha[u_k^n]^\alpha$$

就是使用了一个简写的形式，并不代表$\langle u_k^m\rvert u_q^n\rangle\neq\delta_{m,n}\delta_{k,q}$，其实它们满足的是$\langle u_k^m\rvert u_k^n\rangle=\delta_{m,n}$。到此就先得到了位置算符在占据态上的投影后的结果，因为我们在数值计算的时候，能处理的都是矩阵，basis实际上只是在考虑问题的时候涉及的，所以此时得到的也是一个矩阵

$$[G_k]^{mn}=\langle u_{k+\Delta}^m\rvert u_k^n\rangle$$

可以看到它是哈密顿量不同$k$点处本征态的交叠矩阵，是非幺正的，不过[Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115) 这篇文章中证明了，在热力学极限下面，这个矩阵仍然是幺正的。不过在真实的计算中肯定不能使得系统趋向于热力学极限，为了让其在有限大小$N$的情况下可以是幺正的，那么矩阵$G$进行奇异值分解

$$G=UDV^\dagger$$

其实这个过程就是通常我们处理矩阵对角化的过程，只不过对于一般的矩阵，就先叫做奇异值分解。此时$D$是一个对角矩阵，矩阵$G$是非幺正的特征就是$D$矩阵的对角元素对应的奇异值都是小于1的，所以可以在每一个$k$点定义

$$F=UV^\dagger$$

它是幺正矩阵，将$F$成为在$k$点处的Wilson line element。根据文章所说，在$N\rightarrow\infty$时，$[F_k]^{mn}=[G_k]^{mn}$，这里我就暂时接受这件事情了。为了将投影位置算符$P^{\mathrm{occ}} \hat{x} P^{\mathrm{occ}}$对角化，此时可以得到其对应的本征值问题为

$$P^{\mathrm{occ}} \hat{x} P^{\mathrm{occ}}\rvert\Psi^j\rangle=E^j\rvert\psi^j\rangle$$

根据投影位置算符在前面的表达式，在基矢$\gamma_{n,k}\rvert 0\rangle$下面，上面的本征方程形式为

$$\begin{equation}
\left(\begin{array}{ccccc}
0 & 0 & 0 & \ldots & F_{k_{N}} \\
F_{k_{1}} & 0 & 0 & \ldots & 0 \\
0 & F_{k_{2}} & 0 & \ldots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \ldots & 0
\end{array}\right)\left(\begin{array}{c}
v_{k_{1}} \\
v_{k_{2}} \\
v_{k_{3}} \\
\vdots \\
v_{k_{N}}
\end{array}\right)^{j}=E^{j}\left(\begin{array}{c}
v_{k_{1}} \\
v_{k_{2}} \\
v_{k_{3}} \\
\vdots \\
v_{k_{N}}
\end{array}\right)^{j},
\end{equation}$$

这里离散的$k$点取值为$k_1=0,k_2=\Delta_k,\cdots,k_N=\Delta_k(N-1)$，指标$j\in 1,\cdots N_{\rm occ}$表示占据态。这里先将上面的方程重新写一下

$$
F\mathbf{v}=E\mathbf{v}
$$
这里的

$$
F=\left(\begin{array}{ccccc}
0 & 0 & 0 & \ldots & F_{k_{N}} \\
F_{k_{1}} & 0 & 0 & \ldots & 0 \\
0 & F_{k_{2}} & 0 & \ldots & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
0 & 0 & 0 & \ldots & 0
\end{array}\right)
$$

可以将$F$重复的作用在方程左侧可以得到

$$
F^{N}\mathbf{v}=E\mathbf{v}
$$

这里$N$其实就是这个矩阵的维数，其实自己可以验证一下，对于$F$这种具有特殊形式的矩阵，将其多次作用在一个列向量上之后有

$$
\mathcal{W}_{k_f\leftarrow k_i}\rvert v^j_{k_i}\rangle=(E^j)^{(k_f-k_i)/\Delta_k}\rvert v^j_{k_f}\rangle
$$

这里还是使用了Dirac符号表示了$[v_{k_l}^j]^n$这个列矢量，这里$n\in 1\cdots N_{\rm occ}$，此时可以定义处离散的Wilson line

$$\mathcal{W}_{k_f\leftarrow k_i}=F_{k_f-\Delta_k}F_{k_f-2\Delta_k}\cdots F_{k_f+\Delta_k}F_{k_i}$$

对于足够大的Wilson loop，总是可以将其穿过整个BZ的，此时就可以上面的公式重新表示，令其能覆盖整个BZ

$$\mathcal{W}_{k+2\pi\leftarrow k}\rvert v_k^j\rangle=(E^j)^N\rvert v_k^j\rangle$$

这里的$k$就是计算Wilson loop的起点(base point)，虽然选择不同的$k$可能使得Wilson loop的本征态会有些不同，但是其本征值并不会随着$k$的选择不同而发生变化。而且此时因为Wilson loop是幺正的，那么其本征值可以简单的表示为位相

$$(E^j)^N=e^{i2\pi v^j}$$

即Wilson loop对应的本征方程为

$$\mathcal{W}_{k+2\pi\leftarrow k}\rvert v_k^j\rangle=e^{i2\pi v^j}\rvert v_k^j\rangle$$

这里$j=1\cdots N_{\rm occ}$还是遍历所有的占据态，而$v^j$就是占据态能带对应的Wannier center，因此有多少个占据态就会对应多少个Wannier center，而原胞中体态的极化通常是对所有的占据态进行的，即

$$p=\sum_jv^j$$


将其用Wilson loop表示出来即

	$$P=-\frac{i}{2\pi}\log\det[\mathcal{W}_{k+2\pi\leftarrow k}]$$

	此时来考虑前面提到的矩阵$[G_k]^{mn}=\langle u^m_{k+\Delta_k}\rvert u^n_k\rangle$，这里$\Delta_k=(k_f-k_i)/N$，在热力学极限下$N\rightarrow\infty$，可以将本征态进行展开

	$$\langle u^m_{k+\Delta_k}\rvert=\langle u^m_k\rvert+\Delta_k\partial_k\langle u^m_k\rvert+\cdots$$

	此时可以将矩阵中元素表示为

	$$[G_k]^{mn}=\langle u_k^m\rvert u_k^n\rangle+\Delta_k\langle \partial_k u_k^m\rvert u_k^n\rangle+\cdots$$

利用本征态之间的正交关系

$$\langle u_k^m\rvert u_k^n\rangle=\delta^{mn}\qquad \langle\partial_k u^m_k\rvert u_k^n\rangle=-\langle u^m_k\rvert \partial_k u_k^n\rangle$$

此时可以将其重新表示

$$[G_k]^{mn}=\delta^{mn}-\Delta_k\langle u^m_k\rvert\partial_k u_k^n\rangle=\delta^{mn}-i\Delta_k[\mathcal{A}_k]^{mn}$$

这里Berry联络定义为

$$[\mathcal{A}_k]^{mn}=-i\langle u^m_k\rvert\partial_k\rvert u_k^n\rangle$$

它是一个纯实数量，因此在热力学极限下，可以将Wilson loop表示为

$$\mathcal{W}_{k_f\leftarrow k_i}=\lim_{N\rightarrow\infty}\Pi_{n=1}^N[I-i\Delta_k\mathcal{A}_{k+n\Delta_k}]=\exp[-i\int_{k_i}^{k_f}\mathcal{A}_kdk]$$

此时极化就可以重新表示为

$$p=-\frac{i}{2\pi}\log\det[e^{-i\int_k^{k+2\pi}\mathcal{A}_kdk}]=-\frac{1}{2\pi}\int_k^{k+2\pi}\text{Tr}[\mathcal{A}_k]dk\quad\text{mod}\quad 1$$


所以电子的极化性质就是有占据态空间的Berry位相决定的。


