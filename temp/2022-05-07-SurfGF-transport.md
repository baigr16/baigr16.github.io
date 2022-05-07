---
title: 迭代方法自能计算
tags: Julia transport 
layout: article
license: true
toc: true
key: a20220507a
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
最近在学习输运方面的内容，这篇Blog主要就是整理如何通过迭代格林函数的方式来计算电极与导体耦合之后的自能。不得不说知乎是个学习的好地方，我这里的主要内容也都是在知乎上学习的，自己整理一下再加上一些自己对问题的理解。
{:.info}
<!--more-->
# 自能计算
这里考虑左右两个电极，中间是导体区域

体系的电导$G=\frac{2e^2}{h}T_{LR}$，这里的$T_{LR}$是透射率

$$
T_{LR}=\text{Tr}[\Gamma_LG_{1M}\Gamma_RG_{1M}^\dagger]
$$

公式的来源可以参考Data的介观输运那本书。这里的$\Gamma_L,\Gamma_R$分别是左右电极的线宽，可以通过左右电极的自能
来的到。首先考虑一个由左右两个电极重甲是导体的体系，对应的哈密顿量可以表示为
$$
\left(\begin{array}{ccc}
H_{L} & H_{L C} & 0 \\
H_{C L} & H_{C} & H_{C R} \\
0 & H_{R C} & H_{R}
\end{array}\right)
$$

这里的$H_{LC},H_{RC}$表示的就是导体和电极之间的耦合。那么有了哈密顿量之后，我们就可以来定义格林函数

$$
\left(\begin{array}{ccc}
E I-H_{L} & -H_{L C} & 0 \\
-H_{C L} & E I-H_{C} & -H_{C R} \\
0 & -H_{R C} & E I-H_{R}
\end{array}\right)\left(\begin{array}{ccc}
G_{L L} & G_{L C} & G_{L R} \\
G_{C L} & G_{C C} & G_{C R} \\
G_{R L} & G_{R C} & G_{R R}
\end{array}\right)=\left(\begin{array}{ccc}
I & 0 & 0 \\
0 & I & 0 \\
0 & 0 & I
\end{array}\right)
$$

通常我们需要求解的只是导体的电导$G_{CC}$，所以将格林函数矩阵中包含$G_{CC}$的这一列单独表示出来

$$
\left(\begin{array}{ccc}
E I-H_{L} & -H_{L C} & 0 \\
-H_{C L} & E I-H_{C} & -H_{C R} \\
0 & -H_{R C} & E I-H_{R}
\end{array}\right)\left(\begin{array}{l}
G_{L C} \\
G_{C C} \\
G_{R C}
\end{array}\right)=\left(\begin{array}{l}
0 \\
I \\
0
\end{array}\right)
$$

那么就可以得到$G_{CC}$的表达式

$$
G_{CC}=(EI-H_C-\Sigma_L-\Sigma_R)^{-1}
$$

这里

$$
\Sigma_L=H_{CL}(EI - H_L)^{-1}H_{LC}
$$

表示左电极对应的自能，同样右电极对应的自能为

$$
\Sigma_R=H_{CR}(EI - H_R)^{-1}H_{RC}
$$

之前在学习的时候一直有个疑问，我们要求解的很多物理量最后都与自能有关，但是我们又不知道自能到底是多少，这里就来详细看一下到底要如何计算自能。
{:.warning}

当我们考虑周期系统的时候，原胞内和原胞间的hopping矩阵分别为$H_{00},H_{11}$，此时电极对应的哈密顿量为

$$
H_{R}=\left(\begin{array}{cccc}
H_{00} & H_{01} & 0 & \cdots \\
H_{10} & H_{00} & H_{01} & \cdots \\
0 & H_{10} & H_{00} & \\
\vdots & \vdots & \vdots & \ddots
\end{array}\right)
$$

同样可以对电极定义出格林函数

$$
\left(\begin{array}{cccc}
E I-H_{00} & -H_{01} & 0 & \cdots \\
-H_{10} & E I-H_{00} & -H_{01} & \cdots \\
0 & -H_{10} & E I-H_{00} & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{array}\right)\left(\begin{array}{cccc}
G_{00} & G_{01} & G_{02} & \cdots \\
G_{10} & G_{11} & G_{12} & \cdots \\
G_{20} & G_{21} & G_{22} & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{array}\right)=I
$$

因为我们这里考虑的是个无限大的系统，所以这个矩阵的维度就是无限大，看起来好像我们没有办法去求它对应的本征值和本征态了。但是如果熟悉紧束缚近似，并且通过观察这个矩阵的形式也可以发现
它是具有规律性的，因为这里我们关心的只是最近邻格点之间的hopping，所以这个哈密顿量的形式看起来就是三对角的形式，这也就是它的规律了。因为我们想要求解自能
$\Sigma_R=H_{CR}(EI - H_R)^{-1}H_{RC}$，这里$H_{CR}$表示的是电极和导体之间的耦合，我们在实际考虑的时候也只有在电极与导体接触的那个位置，这个耦合才不为零，所以它的形式为

$$
H_{C R}=\left(\begin{array}{cccc}
H_{01} & 0 & 0 & \ldots \\
0 & 0 & 0 & \ldots \\
0 & 0 & 0 & \ldots \\
\vdots & \vdots & \vdots & \ddots
\end{array}\right)
$$

所以上面的这个耦合矩阵中，只有一个小块才是非零的，而前面给出的格林函数的计算公式中，我们就只需要计算它的第一个块矩阵就可以了(这一个块中包含了导体与电极之间的耦合，其它块中不包含耦合，所以格林函数是可以直接得到的)。这里我们要求姐的格林函数为$G_{00}$，所以将导线格林函数的第一列单独提取出来

$$
\left(\begin{array}{cccc}
E I-H_{00} & -H_{01} & 0 & \cdots \\
-H_{10} & E I-H_{00} & -H_{01} & \cdots \\
0 & -H_{10} & E I-H_{00} & \cdots \\
\vdots & \vdots & \vdots & \ddots
\end{array}\right)\left(\begin{array}{c}
G_{00} \\
G_{10} \\
G_{20} \\
\vdots
\end{array}\right)=\left(\begin{array}{c}
I \\
0 \\
0 \\
\vdots
\end{array}\right)\label{p1}
$$

这又是一个三对角的形式，看到这种形式其实就很愉快，因为它表示的其实就是一个递归关系，但是要从$n\geq2$开始，将这个矩阵表示展开可以得到

$$
-H_{10}G_{n,0}+(EI-H_{00})G_{n + 1,0}-H_{01}G_{n + 2,0}=0
$$

将它表示为矩阵形式

$$
\left(\begin{array}{c}
G_{n+2,0} \\
G_{n+1,0}
\end{array}\right)=\left(\begin{array}{cc}
H_{01}^{-1}\left(E I-H_{00}\right) & H_{01}^{-1} H_{10} \\
I & 0
\end{array}\right)\left(\begin{array}{c}
G_{n+1,0} \\
G_{n, 0}
\end{array}\right)
$$

这不就是量子力学中层间碰到过的转移矩阵的问题了吗，这里将$T$表示为

$$
T=\left(\begin{array}{cc}
H_{01}^{-1}\left(E I-H_{00}\right) & H_{01}^{-1} H_{10} \\
I & 0
\end{array}\right)
$$

可以发现$T$的作用就是实现

$$
\left(\begin{array}{c}
G_{n+1,0} \\
G_{n, 0}
\end{array}\right)\rightarrow\left(\begin{array}{c}
G_{n+2,0} \\
G_{n+1,0}
\end{array}\right)
$$

这样的一个过程，加入将左端从$n=0$开始，右端是$n=N$，那么就是将转移矩阵$T$作用在$n=0$的位置上$N$次。但是乍一看好像还会没有解决问题，因为转移矩阵$T$乘$N$次之后是多少。这个时候物理的思维就显得很重要，我们想要将一个初态通过转移矩阵的方式连续的转移$N$次，要知道每一次转移都是有一定的概率的，而且这个概率肯定是$\leq 1$的，所以当作用$N$次之后，假设$N\rightarrow\infty$，那么转移的几率一定会$p\rightarrow0$。这样就是说转移矩阵满足

$$
T^{N\rightarrow\infty}=0
$$

将方程(\ref{p1})中的第一行乘过去可以得到

$$
(EI - H_{00}-H_{01}G_{10})=I\label{p2}
$$

我们在这里可以假设

$$
G_{10}=\gamma_1G_{00}
$$

将它回代到方程(\ref{p2})中就可以得到$G_{00}$了

$$
G_{00}=(EI-H_{00}-H_{01}\gamma_1)
$$

所以这个时候需要解决的问题就是怎么求解$\gamma_1$。

将方程(\ref{p1})中的第二行乘进去，并且此时假设$G_{20}=\gamma_2G_{00}$，通过代数运算可以得到

$$
\gamma_1=g_1H_{10}+g_1H_{01}\gamma_2\quad g_1=(EI - H_{00})^{-1}
$$

好了，$\gamma_1$求解出来了，但是它又包含了$\gamma_2$，那么我们就需要再求解$\gamma_2$。同样的过程可以对公式(\ref{p1})中后面的每一行都进行，可以得到

$$
\begin{array}{c}
\gamma_{2}=g_{2} H_{10} g_{1} H_{10}+g_{2} H_{01} \gamma_{3} \\
\gamma_{3}=g_{3} H_{10} g_{2} H_{10} g_{1} H_{10}+g_{3} H_{01} \gamma_{4} \\
\gamma_{4}=g_{4} H_{01} g_{3} H_{10} g_{2} H_{10} g_{1} H_{10}+g_{4} H_{01} \gamma_{5} \\
\ldots
\end{array}\label{p3}
$$

这里的$g_i$的形式会比较复杂，可以自己推导试试，我也是在看过知乎上的推导实际操作了一番，发现形式上还是有点繁琐，但是没关系，我们需要的只是一个规律而已，通过推导就可以整理出上面的表达式。但是又发现问题了，这个表达式没有尽头求解$\gamma_i$会耦合出来$\gamma_{i+1}$，那么此时格林函数就显的很智慧了。因为格林函数的本质就是点与点之间的关联，我们自然可以认为
更远位置处的关联是零，那么就可以公式(\ref{p3})在某一阶直接截断，即零$\gamma_i=0$，那么自然问题就解耦了。这里就抄一下知乎上的结果，将$\gamma_5=0$可以得到

$$
\begin{array}{l}
\gamma_{4}=g_{4} H_{01} g_{3} H_{10} g_{2} H_{10} g_{1} H_{10} \\
\gamma_{3}=g_{3} H_{10} g_{2} H_{10} g_{1} H_{10}+g_{3} H_{01} g_{4} H_{01} g_{3} H_{10} g_{2} H_{10} g_{1} H_{10} \\
\gamma_{2}=g_{2} H_{10} g_{1} H_{10}+g_{2} H_{01} g_{3} H_{10} g_{2} H_{10} g_{1} H_{10}+g_{2} H_{01} g_{3} H_{01} g_{4} H_{01} g_{3} H_{10} g_{2} H_{10} g_{1} H_{10} \\
\gamma_{1}=g_{1} H_{10}+g_{1} H_{01} g_{2} H_{10} g_{1} H_{10}+g_{1} H_{01} g_{2} H_{01} g_{3} H_{10} g_{2} H_{10} g_{1} H_{10}+g_{1} H_{01} g_{2} H_{01} g_{3} H_{01} g_{4} H_{01} g_{3} H_{10} g_{2} H_{10} g_{1} H_{10}
\end{array}
$$

好长的式子，但是没关系，它是有规律可循的，我们可以通过迭代的方式来计算，代码如下

```julia
# 通过迭代的方式求解输运中的表面格林函数
#------------------------------------------
using LinearAlgebra,DelimitedFiles,PyPlot
#-----------------------------------------
function matset()
    e0::Float64 = 4
    t0::Float64 = 1.0
    hn::Int64 = 4
    h00 = zeros(ComplexF64,hn,hn)
    h01 = zeros(ComplexF64,hn,hn)
    h00[1,1] = e0
    h00[1,2] = t0
    h00[2,2] = e0
    h00[2,1] = t0
    h00[3,3] = e0
    h00[3,4] = t0
    h00[4,4] = e0
    h00[4,3] = t0
    #---------------
    h01[1,1] = t0
    h01[2,2] = t0
    h01[3,3] = t0
    h01[4,4] = t0
    #---------------
    return h00,h01
end
#-----------------------------------------------
function SurfGF()
    err::Float64 = 1.0 # 控制误差
    ha,a0 = matset()
    H00,H01 = ha,a0
    hb,b0 = ha,a0
    hn = size(ha)[1] # 获取矩阵维度
    eye = zeros(ComplexF64,hn,hn)
    g0 = zeros(ComplexF64,hn,hn)
    eng = 1.0 + 0.01im # 给定能量计算格林函数,这里需要加入小虚部,否则函数会出现发散的情况
    for i0 in 1:hn
        eye[i0,i0] = 1 #构造单位矩阵
    end
    #-----------
    cont = 0
    while err > 10^(-9)
    # while cont < 20
        h0 = g0
        g0 = inv(eng * eye - hb)
        ha += (a0 * g0 * b0)
        hb += (a0 * g0 * b0 + b0 * g0 * a0)
        a0 = a0 * g0 * a0
        b0 = b0 * g0 * b0
        cont += 1
        err = abs(sum(g0 - h0))
    end
    g1 = inv(eng * eye - ha) #计算得到的就是表面格林函数
    g2 = H01 * g1 * H01' # 自能
    g3 = im * (g2 - conj(g2)) # 展宽函数
    println(g2)
    return g1,g2,g3
end
#-------------------------------------------------
@time g1,g2,g3 = SurfGF()
```
通过计算发现这个计算还是非常快的，这也体现了格林函数方法的高效性。

# 参考
- 1.[输运问题总结与练习(一):二端体系](https://zhuanlan.zhihu.com/p/269595149)