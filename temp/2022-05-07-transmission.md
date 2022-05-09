---
title: 格林函数方法计算透射率
tags: Julia transport 
layout: article
license: true
toc: true
key: a20220507b
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
前面的一篇Blog主要是整理了一下怎么通过迭代的方式来得到电极对应的自能。这里就再来看一下导体区域的格林函数怎么得到，并通过与前面自能的结合来计算透射率。
{:.info}
<!--more-->
# 前言

通过上一篇Blog已经可以计算得到电极表面格林函数了$g_r$，那么相应的就可以得到自能

$$
\Sigma_r=H_{01}g_rH_{10}
$$

对应的展宽函数也可以得到

$$
\Gamma_r=i*(\Sigma_r-\Sigma_r^\dagger)
$$

透射率的计算公式为

$$
T_{LR}=\text{Tr}[\Gamma_LG_C\Gamma_RG_C^\dagger]
$$

我们已经知道展宽函数其实只有一个块矩阵是非零的，所以我们这里先对$\Gamma_LG_C\Gamma_RG_C^\dagger$这个矩阵进行分析，将它展开可以得到

$$
\begin{array}{l}
\Gamma_{L} G_{C} \Gamma_{R} G_{C}^{\dagger}=\left(\begin{array}{cccc}
\Gamma_{L} & 0 & \cdots & 0 \\
0 & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 0
\end{array}\right)\left(\begin{array}{cccc}
G_{11} & G_{12} & \cdots & G_{1 M} \\
G_{21} & G_{22} & \cdots & G_{2 M} \\
\vdots & \vdots & \ddots & \vdots \\
G_{M 1} & G_{M 2} & \cdots & G_{M M}
\end{array}\right)\left(\begin{array}{cccc}
0 & 0 & \cdots & 0 \\
0 & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & \Gamma_{R}
\end{array}\right)\left(\begin{array}{ccccc}
G_{11} & G_{12} & \cdots & G_{1 M} \\
G_{21} & G_{22} & \cdots & G_{2 M} \\
\vdots & \vdots & \ddots & \vdots \\
G_{M 1} & G_{M 2} & \cdots & G_{M M}
\end{array}\right)\\
=\left(\begin{array}{cccc}
\Gamma_{L} G_{11} & \Gamma_{L} G_{12} & \cdots & \Gamma_{L} G_{1 M} \\
0 & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 0
\end{array}\right)\left(\begin{array}{cccc}
0 & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & 0 \\
\Gamma_{R} G_{1 M}^{\dagger} & \Gamma_{R} G_{2 M}^{\dagger} & \cdots & \Gamma_{R} G_{M M}^{\dagger}
\end{array}\right)\\
=\left(\begin{array}{cccc}
\Gamma_{L} G_{1 M} \Gamma_{R} G_{1 M}^{\dagger} & 0 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
\Gamma_{L} G_{1 M} \Gamma_{R} G_{M-1, M}^{\dagger} & 0 & \cdots & 0 \\
\Gamma_{L} G_{1 M} \Gamma_{R} G_{M M}^{\dagger} & 0 & \cdots & 0
\end{array}\right)
\end{array}
$$

可以发现因为之前从电极得到的展宽函数只有一个小块是非零的，此时我们得到的结果就只有一列是非零的，那么$T_{LR}=\text{Tr}[\Gamma_LG_C\Gamma_RG_C^\dagger]$就变成了只对结果中的第一个块矩阵求trace即可，可就是说现在唯一不知道的就是中心散射区域的格林函数$G_{1M}$，下面就来求解它。首先还是先把中心散射区的格林函数表示出来

$$
\left(\begin{array}{cccc}
E I-H_{00}-\Sigma_{L} & -H_{01} & \cdots & 0 \\
-H_{10} & E I-H_{00} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & E I-H_{00}-\Sigma_{R}
\end{array}\right)\left(\begin{array}{cccc}
G_{11} & G_{12} & \cdots & G_{1 M} \\
G_{21} & G_{22} & \cdots & G_{2 M} \\
\vdots & \vdots & \ddots & \vdots \\
G_{M 1} & G_{M 2} & \cdots & G_{M M}
\end{array}\right)=I
$$

因为我们只想得到$G_{1M}$，所以就将最后一列格林函数单独拿出来

$$
\left(\begin{array}{cccc}
E I-H_{00}-\Sigma_{L} & -H_{01} & \cdots & 0 \\
-H_{10} & E I-H_{00} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & E I-H_{00}-\Sigma_{R}
\end{array}\right)\left(\begin{array}{c}
G_{1 M} \\
G_{2 M} \\
\vdots \\
G_{M M}
\end{array}\right)=\left(\begin{array}{c}
0 \\
0 \\
\vdots \\
I
\end{array}\right)\label{p1}
$$

这个形式又很熟悉，三对角的形式。其实线性代数的就是方程组的问题，我们就来通过消元法来对矩阵中的元素进行操作。因为方程最右端只有最后一个元素是非零的，所以根据行操作来进行消元，我们消元的依据就是就想让矩阵变成一个对角的矩阵。
所以先对第一行乘以$H_{10}G_{1}$再加到第二行上面，这里$G_1=(EI-H_{00}-\Sigma_L)^{-1}$，那么可以得到

$$
\left(\begin{array}{cccc}
G_{1}^{-1} & -H_{01} & \cdots & 0 \\
0 & E I-H_{00}-H_{10} G_{1} H_{01} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & E I-H_{00}-\Sigma_{R}
\end{array}\right)
$$

再来定义$G_2=(EI-H_{00}-H_{10}G_1H_{01})^{-1}$，将第二行乘以$H_{10}G_2$再加到第三行上面可以得到

$$
\left(\begin{array}{ccccccc}
G_{1}^{-1} & -H_{01} & 0 & 0 & \cdots & 0 & 0 \\
0 & G_{2}^{-1} & -H_{01} & 0 & \cdots & 0 & 0 \\
0 & 0 & E I-H_{00}-H_{10} G_{2} H_{01} & -H_{01} & \cdots & 0 & 0 \\
0 & 0 & -H_{10} & E I-H_{00} & \cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & 0 & 0 & \cdots & E I-H_{00} & -H_{01} \\
0 & 0 & 0 & 0 & \cdots & -H_{10} & E I-H_{00}-\Sigma_{R}
\end{array}\right)
$$

将这个过程依次往下进行，最后就可以得到

$$
\left(\begin{array}{ccccccc}
G_{1}^{-1} & -H_{01} & 0 & 0 & \cdots & 0 & 0 \\
0 & G_{2}^{-1} & -H_{01} & 0 & \cdots & 0 & 0 \\
0 & 0 & G_{3}^{-1} & -H_{01} & \cdots & 0 & 0 \\
0 & 0 & 0 & G_{4}^{-1} & \cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & 0 & 0 & \cdots & G_{M-1}^{-1} & -H_{01} \\
0 & 0 & 0 & 0 & 0 & 0 & G_{M}^{-1}
\end{array}\right)\label{p2}
$$

这里有$G_i=(EI-H_{00}-H_{10}G_{i-1}H_{01})^{-1}$，要注意在$i=1$和$i=M$的时候，因为这个位置是有电极存在的，所以此时这个格林函数中还需要考虑上电子带来的自能。我们通过比较方程(\ref{p1})和(\ref{p2})可以发现此时

$$
G_{MM}=G_{M}
$$

但是还不够，我们想要得到的是$G_{1M}$，那么我们再从最后一行往上进行消元，将最后一行乘以$H_{01}G_M$加到倒数第二行上面，因为此时(\ref{p1})最后端的列向量最后一个元素是非零的，所以此时它也会发生变化，要加上同样的数。倒数第二行乘以$G_{M-1}H_{01}$加到倒数第三行上面，此时行向量就变成了

$$
H_{01}G_{M-1}H_{01}G_M
$$

到这里格林函数的结构就变得比较明显了，随着所有的非对角元素通过上面的消元方式被消去，(\ref{p1})左右端的列向量上面的值变为

$$
G_{ij}=H_{01}G_{i+1}H_{01}G_{i+2}H_{01}\cdots G_{M-1}H_{01}G_M
$$

还是要注意，格林函数开始的最左端和结束的最右端是要考虑上自能的，而中间的格林函数则不需要。整理之后可以得到矩阵方程为

$$
\left(\begin{array}{cccc}
G_{1}^{-1} & 0 & \cdots & 0 \\
0 & G_{2}^{-1} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & G_{M}^{-1}
\end{array}\right)\left(\begin{array}{c}
G_{1 M} \\
G_{2 M} \\
\vdots \\
G_{M M}
\end{array}\right)=\left(\begin{array}{c}
H_{01} G_{2} H_{01} G_{3} H_{01} \cdots G_{M-1} H_{01} G_{M} \\
H_{01} G_{3} H_{01} \cdots G_{M-1} H_{01} G_{M} \\
\vdots \\
I
\end{array}\right)
$$

此时中心散射区域的格林函数$G_{1M}$就显而易见了，那么我们就可以取求解对应的透射率了。

# 散射区格林函数
```julia
using LinearAlgebra,PyPlot
#-------------------------------------------
function matset()
    hn::Int64 = 10 # 格点数目
    w::Float64 = 1.0
    v::Float64 = 4.0
    hnn::Int64 = 4 #晶格内部自由度
    eye = Diagonal(repeat([1],outer  = hn)) # 构造单位矩阵
    d1 = Bidiagonal(repeat([0],outer  = hn),repeat([1],outer  = hn - 1),:U) # 负对角线
    d2 = Bidiagonal(repeat([0],outer  = hn),repeat([1],outer  = hn - 1),:L) # 负对角线
    h00 = zeros(Float64,hnn,hnn)
    tx = zeros(Float64,hnn,hnn)
    ty = zeros(Float64,hnn,hnn)
    #------------------------------
    h00[1,2] = w
    h00[1,3] = w
    h00[2,4] = w
    h00[2,1] = w
    h00[3,1] = w
    h00[3,4] = w
    h00[4,2] = w
    h00[4,3] = w
    #--------------------
    tx[1,3] = v
    tx[2,4] = v
    ty[1,2] = v
    ty[3,4] = v
    #--------------------
    H00 = kron(eye,h00) + kron(d1,tx) + kron(d2,tx')
    H01 = kron(eye,ty)
    num = 50 #k点撒点数目
    k = range(-1,1,length = num)
    eng = zeros(Float64,hn * hnn,num)
    for i0 in 1:num
        ham = H00 + H01 * exp(1im * k[i0] * pi) + H01' * exp(-1im * k[i0] * pi)
        eng[:,i0] = map(real,eigvals(ham))
    end
    figure(figsize=(10,8))
    for i0 in 1:hn * hnn
        plot(k,eng[i0,:])
    end
    # println(eng[:,1])
    savefig("test.png",bbox_inches="tight",dpi=300)  # 保存作图文件
    return k,eng
end
#-----------------------------------------
@time matset()
```
# 电导率计算
```julia
using LinearAlgebra,PyPlot
#-------------------------------------------
function matset1()
    hn::Int64 = 10 # 格点数目
    w::Float64 = 1.0
    v::Float64 = 1.0
    hnn::Int64 = 4 #晶格内部自由度
    eye = Diagonal(repeat([1],outer  = hn)) # 构造单位矩阵
    d1 = Bidiagonal(repeat([0],outer  = hn),repeat([1],outer  = hn - 1),:U) # 负对角线
    d2 = Bidiagonal(repeat([0],outer  = hn),repeat([1],outer  = hn - 1),:L) # 负对角线
    h00 = zeros(Float64,hnn,hnn)
    tx = zeros(Float64,hnn,hnn)
    ty = zeros(Float64,hnn,hnn)
    #------------------------------
    h00[1,2] = w
    h00[1,3] = w
    h00[2,1] = w
    h00[2,4] = w
    h00[3,1] = w
    h00[3,4] = w
    h00[4,2] = w
    h00[4,3] = w
    #--------------------
    tx[1,3] = v
    tx[2,4] = v
    ty[1,2] = v
    ty[3,4] = v
    #--------------------
    H00 = kron(eye,h00) + kron(d1,tx) + kron(d2,tx')
    H01 = kron(eye,ty)
    return H00,H01
end
#-------------------------------------------
function matset()
    hn::Int64 = 10 # 格点数目
    u::Float64 = -1.5
    hnn::Int64 = 2 #晶格内部自由度
    eye = Diagonal(repeat([1],outer  = hn)) # 构造单位矩阵
    d1 = Bidiagonal(repeat([0],outer  = hn),repeat([1],outer  = hn - 1),:U) # 负对角线
    d2 = Bidiagonal(repeat([0],outer  = hn),repeat([1],outer  = hn - 1),:L) # 负对角线
    h00 = zeros(ComplexF64,hnn,hnn)
    tx = zeros(ComplexF64,hnn,hnn)
    ty = zeros(ComplexF64,hnn,hnn)
    #------------------------------
    h00[1,1] = u
    h00[2,2] = -u
    #--------------------
    tx[1,1] = 1/2
    tx[2,2] = -1/2
    tx[1,2] = 1im/2
    tx[2,1] = 1im/2

    ty[1,1] = 1/2
    ty[2,2] = -1/2
    ty[1,2] = 1/2
    ty[2,1] = -1/2
    #--------------------
    H00 = kron(eye,h00) + kron(d1,tx) + kron(d2,tx')
    H01 = kron(eye,ty)
    return H00,H01
end
#-------------------------------------------------------------------
function self(H00,H01,ef)
# 给出hopping矩阵计算对应的自能,这里需要输入hopping以及onsite矩阵
    hn = size(H00)[1] # 得到哈密顿量矩阵的维度
    err = 2.0 #统计误差
    eye = Diagonal(repeat([1],outer  = hn)) # 构造单位矩阵
    g0 = zeros(ComplexF64,hn,hn)
    ha,a0 = H00,H01 # 得到hopping矩阵
    hb,b0 = ha,a0'
    cont = 1
    while err > 10^(-6)
    # while cont < 20
        h0 = g0
        g0 = inv(ef * eye - hb)
        ha += (a0 * g0 * b0)
        hb += (a0 * g0 * b0 + b0 * g0 * a0)
        a0 = a0 * g0 * a0
        b0 = b0 * g0 * b0
        cont += 1
        err = abs(sum(g0 - h0))
    end
    g1 = inv(ef * eye - ha) #计算得到的就是表面格林函数
    g2 = H01 * g1 * H01' # 自能
    return g2  # 最后返回的就是自能函数
end
#-------------------------------------------------------------------
function conGF()
    # 计算导体区域的格林函数
    H00,H01 = matset()
    H10 = H01'
    hn::Int64 = size(H00)[1] #矩阵维度
    num = 100
    englist = range(-2,2,length = num)
    trans = []
    eye = Diagonal(repeat([1],outer  = hn)) # 构造单位矩阵
    for i0 in 1:num
        eng = englist[i0] + 0.00001*1im
        selfR = self(H00,H01,eng) #左端自能
        selfL = self(H00,H10,eng) #右端自能
        #--------------------------------------
        GRii = inv(eng * eye - H00 - selfL) #最左端的需要考虑到左端电极的自能
        GR1j = GRii
        # 迭代方式计算中心区域GF
        for k in 2:hn #体系的另外一个维度的小，这里就先设置成正方的
            GRii = inv(eng * eye - H00 - H10 * GRii * H01)
            GR1j = GR1j * H01 * GRii
        end
        # 最右端的要考虑右端电极自能
        GRii = inv(inv(GRii) - selfR)
        # GR1j = GR1j +  GR1j * selfR * GRii

        TR = 1im * (selfR - selfR') #左端线宽函数
        TL = 1im * (selfL - selfL') # 右端线宽

        append!(trans,real(tr(TL * GR1j * TR * GR1j')))
    end
    plot(englist,trans)
    savefig("a.png",bbox_inches="tight",dpi=300)  # 保存作图文件
    return englist,trans

end
#-------------------------------------------------------------------
function main()
    a1,a2 = conGF()
end
#--------------------------------------------------------------------
@time main()
```



# 参考
- 1.[输运问题总结与练习(一):二端体系](https://zhuanlan.zhihu.com/p/269595149)