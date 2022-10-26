---
title: 实空间计算电四极矩
tags:   Topolog Julia
layout: article
license: true
toc: true
key: a20221023
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
这里整理一下如何在实空间中计算高阶拓扑相的电四极矩，如果在体系具有平移对称性的时候，那么自然可以在动量空间中进行计算，但是如果系统中存在杂质或者无序，此时平移对称性被破坏了，然而通常我们考虑的都是强拓扑，所以在实空间中仍然是可以计算其对应的拓扑不变量的。
{:.info}
<!--more-->
# 实空间电四极矩
这里暂时先不讨论实空间中如何通过位置算符得到电四极矩，就直接给出计算表达式

$$q_{xy}=\frac{1}{2\pi}\text{Im}\log[\det(U^\dagger\hat{Q}U)\sqrt{\det(\hat{Q}^\dagger)}]$$

公式中

$$\hat{Q}\equiv \exp(i2\pi\hat{q}_{xy})\qquad \hat{q}_{xy}\equiv\hat{x}\hat{y}/(L_xL_y)$$

这里$\hat{x}(\hat{y})$是沿着$x(y)$方向的位置算符，$L_{x,y}$则是格点尺寸，且

$$\sqrt{\det(\hat{Q}^\dagger)}=\exp[-i\pi\text{Tr}\hat{q}_{xy}]$$

这里$U$矩阵由占据态本征态构成，$U^\dagger U$即为占据态空间的投影算符。考虑[BBH模型](https://link.aps.org/doi/10.1103/PhysRevLett.125.166801)

$$\begin{aligned}
H_{q}(\mathbf{k})=& t \sin k_{y} \gamma_{1}+\left[t_{y}+t \cos k_{y}\right] \gamma_{2} \\
&+t \sin k_{x} \gamma_{3}+\left[t_{x}+t \cos k_{x}\right] \gamma_{4},\label{q1}
\end{aligned}$$

# 分析
这里整理一下自己在写程序时候的思路，因为是要考虑实空间中的哈密顿量，所以首先就需要对动量空间中的哈密顿量\eqref{q1}进行Fourier变换，得到实空间中的哈密顿量，即

$$C^\dagger_{\mathbf{k}}=\frac{1}{\sqrt{N}}\sum_{\mathbf{R}}C^\dagger_{\mathbf{R}}e^{i\mathbf{k}\cdot\mathbf{R}}$$

在变换到实空间的时候，通常会选择两种不同的基矢

$$\Psi^\dagger_1=(C^\dagger_{1a\uparrow},C^\dagger_{1a\downarrow},C^\dagger_{1b\uparrow},C^\dagger_{1b\downarrow};C^\dagger_{2a\uparrow},C^\dagger_{2a\downarrow},C^\dagger_{2b\uparrow},C^\dagger_{2b\downarrow};\cdots,C^\dagger_{Na\uparrow},C^\dagger_{Na\downarrow},C^\dagger_{Nb\uparrow},C^\dagger_{Nb\downarrow})$$

$$\begin{align}\Psi_2^\dagger=(&C^\dagger_{1a\uparrow},C^\dagger_{2a\uparrow},\cdots,C^\dagger_{Na\uparrow};C^\dagger_{1a\downarrow},C^\dagger_{2a\downarrow},\cdots,C^\dagger_{Na\downarrow};\\ &C^\dagger_{1b\uparrow},C^\dagger_{2b\uparrow},\cdots,C^\dagger_{Nb\uparrow};C^\dagger_{1b\downarrow},C^\dagger_{2b\downarrow},\cdots,C^\dagger_{Nb\downarrow})\end{align}$$

可以看到这两种基矢的选择，对应的哈密顿量其实就会不同，虽然它们之间相差的只是一个幺正变换，但是此时需要注意的是，不同的基矢选择，在构建位置算符的时候，也要对应的进行修改。

假如我现在选择第一种basis

$$\Psi^\dagger_1=(C^\dagger_{1a\uparrow},C^\dagger_{1a\downarrow},C^\dagger_{1b\uparrow},C^\dagger_{1b\downarrow};C^\dagger_{2a\uparrow},C^\dagger_{2a\downarrow},C^\dagger_{2b\uparrow},C^\dagger_{2b\downarrow};\cdots,C^\dagger_{Na\uparrow},C^\dagger_{Na\downarrow},C^\dagger_{Nb\uparrow},C^\dagger_{Nb\downarrow})$$

此时位置算符$\hat{q}_{xy}$的形式应该为

$$\hat{q}_{xy}=\frac{\hat{x}\hat{y}}{L_xL_y}=(\frac{1}{N},\frac{1}{N},\frac{1}{N},\frac{1}{N};\frac{2}{N},\frac{2}{N},\frac{2}{N},\frac{2}{N};\cdots,\frac{N}{N},\frac{N}{N},\frac{N}{N},\frac{N}{N})$$

这里假设实空间中$x$方向和$y$方向的格点数目分别是$L_x$和$L_y$，那么对应的$N=L_x\times L_y$。

类似的，如果选择第二种基矢

$$\begin{align}\Psi_2^\dagger=(&C^\dagger_{1a\uparrow},C^\dagger_{2a\uparrow},\cdots,C^\dagger_{Na\uparrow};C^\dagger_{1a\downarrow},C^\dagger_{2a\downarrow},\cdots,C^\dagger_{Na\downarrow};\\ &C^\dagger_{1b\uparrow},C^\dagger_{2b\uparrow},\cdots,C^\dagger_{Nb\uparrow};C^\dagger_{1b\downarrow},C^\dagger_{2b\downarrow},\cdots,C^\dagger_{Nb\downarrow})\end{align}$$

此时位置算符$\hat{q}_{xy}$的形式为

$$\begin{align}\hat{q}_{xy}=(&\frac{1}{N},\frac{2}{N},\cdots,\frac{N}{N}; \frac{1}{N},\frac{2}{N},\cdots,\frac{N}{N}\\ &\frac{1}{N},\frac{2}{N},\cdots,\frac{N}{N};\frac{1}{N},\frac{2}{N},\cdots,\frac{N}{N}) \end{align}$$

因此，需要根据自己实空间哈密顿量basis的选择，进一步再构建正确的位置算符。

# 代码实现
我个人在构建实空间哈密顿量的时候，习惯的方式是$\Psi_2^\dagger$这种basis，在二次量子化表示下，体系的哈密顿量表示为

$$\hat{H}(\mathbf{k})=\sum_\mathbf{k}\Psi^\dagger_\mathbf{k}H(\mathbf{k})\Psi_\mathbf{k}$$

在利用程序处理矩阵对角化的时候，其实就是在处理$H(\mathbf{k})$这个矩阵。还是利用[BBH](https://link.aps.org/doi/10.1103/PhysRevB.96.245115)作为例子来进行研究，首先实空间中的哈密顿量构造为

```julia
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles,Printf
# --------------------------------------
@everywhere function boundary(xn::Int64,yn::Int64)
    len2::Int64 = xn*yn
    bry = zeros(Int64,4,len2)
    for iy in 1:yn
        for ix in 1:xn
            i = (iy - 1)*xn + ix
            bry[1,i] = i + 1
            if ix == xn
                bry[1,i] = bry[1,i] - xn
            end
            bry[2,i] = i - 1
            if ix == 1
                bry[2,i] = bry[2,i] + xn
            end
            bry[3,i] = i + xn
            if iy == yn 
                bry[3,i] = bry[3,i] - len2
            end
            bry[4,i] = i - xn
            if iy == 1
                bry[4,i] = bry[4,i] + len2
            end
        end
    end
    return bry
end
#-------------------------------------
@everywhere function pauli()
    s0 = zeros(ComplexF64,2,2)
    s1 = zeros(ComplexF64,2,2)
    s2 = zeros(ComplexF64,2,2)
    s3 = zeros(ComplexF64,2,2)
    #----
    s0[1,1] = 1
    s0[2,2] = 1
    #----
    s1[1,2] = 1
    s1[2,1] = 1
    #----
    s2[1,2] = -im
    s2[2,1] = im
    #-----
    s3[1,1] = 1
    s3[2,2] = -1
    #-----
    return s0,s1,s2,s3
end
#---------------------------------------
@everywhere function gamma()
    s0,sx,sy,sz = pauli()
    g1 = -kron(sy,sx) 
    g2 = -kron(sy,sy) 
    g3 = -kron(sy,sz) 
    g4 = kron(sx,s0) 
    g5 = kron(sz,s0) # 微扰
    return g1,g2,g3,g4,g5
end
#---------------------------------------
@everywhere function matset(xn::Int64,yn::Int64,gamx::Float64)
    # gamx::Float64 = 0.5  
    lamx::Float64 = 0.5  
    gamy::Float64 = gamx
    lamy::Float64 = lamx
    d0::Float64 = 0.0
    hn::Int64 = 4
    N::Int64 = xn*yn*hn
    len2::Int64 = xn*yn
    ham = zeros(ComplexF64,N,N)
    g1,g2,g3,g4,g5 = gamma()
    bry = boundary(xn,yn)
    #---
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            for i1 in 0:hn -1 
                for i2 in 0:hn - 1
                    ham[i0 + len2*i1,i0 + len2*i2] = gamx*g4[i1 + 1,i2 + 1] + gamy*g2[i1 + 1,i2 + 1] + d0*g5[i1 + 1,i2 + 1]
                    if ix != xn
                        ham[i0 + len2*i1,bry[1,i0] + len2*i2] = lamx/2.0*g4[i1 + 1,i2 + 1] + lamx/(2*im)*g3[i1 + 1,i2 + 1]
                    end
                    if ix != 1
                        ham[i0 + len2*i1,bry[2,i0] + len2*i2] = lamx/2.0*g4[i1 + 1,i2 + 1] - lamx/(2*im)*g3[i1 + 1,i2 + 1]
                    end
                    if iy != yn
                        ham[i0 + len2*i1,bry[3,i0] + len2*i2] = lamy/2.0*g2[i1 + 1,i2 + 1] + lamy/(2*im)*g1[i1 + 1,i2 + 1]
                    end
                    if iy != 1
                        ham[i0 + len2*i1,bry[4,i0] + len2*i2] = lamy/2.0*g2[i1 + 1,i2 + 1] - lamy/(2*im)*g1[i1 + 1,i2 + 1]
                    end
                end
            end
        end
    end
    #-----------------------------------------
    if ishermitian(ham)
        val,vec = eigen(ham)
    else
        println("Hamiltonian is not hermitian")
        # break
    end
    # fx1 = "juliaval-" * string(h0) * ".dat"
    fx1 = "eigval.dat"
    f1 = open(fx1,"w")
    ind = (a->(@sprintf "%5.2f" a)).(range(1,length(val),length = length(val)))
    val2 = (a->(@sprintf "%15.8f" a)).(map(real,val))
    # writedlm(f1,map(real,val),"\t")
    writedlm(f1,[ind val2],"\t")
    close(f1)
    return map(real,val),vec
end
```

接下来就需要构建在$\Psi^\dagger_2$基矢下面的位置算符。在利用`Julia`编程的时候，这两种位置算符的basis都是比较容易构造的
```julia
xn::Int64 = 10 # 格点数目
hn::Int64 = 4  # 动量空间哈密顿量维度
qxy = repeat(1:xn*xn,inner = hn)/(xn*xn) # Basis 1
qxy = repeat(1:xn*xn,outer = hn)/(xn*xn) # Basis 2
```