---
title: 边界极化
tags:  Topology
layout: article
license: true
toc: true
key: a20221011
pageview: true
cover: /assets/images/topology/edge-pora.png
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
\rvert \Psi^j_{R_x}\rangle = \frac{1}{\sqrt{N_x}} \sum_{n=1}^{N_{occ} \times N_y}\sum_{k_x} \left[ \nu^j_{k_x} \right]^n e^{-i k_x R_x} \gamma^\dagger_{n,k_x}\rvert 0\rangle,
\label{eq3}
\end{align}$$

这里$j\in 1\cdots N_{\rm occ}\times N_y, R_x\in 1\cdots N_x$, $[v_{k_x}^j]^n$是第$j$个Wilson loop本征态$\rvert v_{k_x}^j\rangle$的第$n$个分量, $\gamma^\dagger_{n,k_x}$如公式\eqref{eq4}所示. 在得到了Hybrid Wannier函数之后, 就可以得到其对应的几率密度

$$\begin{align}
\rho^{j,R_x}(R_y) &= \sum_{R'_x,\alpha} \langle \Psi^j_{R_x}\rvert \phi^{R_y, \alpha}_{R'_x}\rangle \langle\phi^{R_y, \alpha}_{R'_x}\rvert\Psi^j_{R_x}\rangle\nonumber\\
&=\frac{1}{N_x} \sum_{k_x, \alpha} \left| [u^n_{k_x}]^{R_y, \alpha}[\nu^j_{k_x}]^n\right|^2\label{eqq}
\end{align}$$

通过$\rho^{j,R_x}(R_y)$就可以进一步在$y$方向的格点上解析出极化在$x$方向的分量. 公式\eqref{eqq}中的符号含义可能需要搞清楚，其中$[u_{k_x}^n]^{R_y,\alpha}$表示的是沿着$y$方向开边界之后，哈密顿量本征态，$n$表示的是第$n$个占据态，而在$[v_{k_x}^j]^n$这里的$n$表示的则是第$j$个Wilson loop本征态的第$n$个分量，所以这里的符号就有点让人迷茫。程序计算的话就是
```julia
@everywhere function rho(yn::Int64,h0::Float64,dir::Int64)
    # 计算所有格点上的极化分布
    # yn 开边界方向原胞数. h0 层间耦合强度
    hn::Int64 = 8
    N::Int64 = yn*hn
    Nocc::Int64 = Int(N/2)
    kn::Int64 = 50
    re1 = zeros(Float64,yn)
    # re1 = SharedArray(zeros(Float64,yn))
    ham_wan = Wannier(h0,yn,dir) # 得到Wannier哈密顿量
    val2,vec2 = eigen(ham_wan)
    val2 = map(real,val2)/(2*pi)
    for ik = -kn:kn
        kx = ik*pi/kn
        if dir == 0
            ham = openx(h0,yn,kx)
        else
            ham = openy(h0,yn,kx)
        end
        val1,vec1 = eigen(ham)
        for j1 = 1:Nocc # Wannier哈密顿量要对所有本征矢量求和
            for n1 = 1:Nocc # 哈密顿量要对所有占据态求和
                for ibeta = 1:hn # 对一个格点上的本征矢量所有分量求和
                    for iy = 1:yn # 计算每个格点上的边界极化
                        re1[iy] = re1[iy] + abs(vec1[(iy - 1)*hn + ibeta,n1]*vec2[n1,j1])^2*val2[j1]
                    end
                end
            end 
        end 
    end
    return re1/(2*kn + 1)
end
```
特别有用的一点就是可以得到在边界上$R_y=1,N_y$,是否存在具有的Hybrid Wannier函数. 通过几率密度$\rho^{j,R_x}(R_y)$可以计算极化在$x$方向的分量

$$\begin{align}
p_x(R_y) = \sum_j \rho^j(R_y) \nu_x^j 
\label{eq5}
\end{align}$$

从而就可以得到极化在$y$边界上是如何分布的, 通过变化格点$R_y$即可. 为了得到边界极化, 可以从边界中点到边界的上对极化分布进行求和

$$p_x^{\rm edge, y}=\sum_{i=1}^{N_y/2}p_x(i_y)$$

从而就可以得到边界极化. 上面的这些过程同样也可以是沿着$y$方向开边界进行计算的。

![png](/assets/images/topology/edge-pora.png)

# 参考
- [Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)