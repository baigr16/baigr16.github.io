---
title: BBH模型实空间计算电四极矩
tags: Topology Julia
layout: article
license: true
toc: true
key: a20220525
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
这篇Blog整理一下怎么在实空间计算电四极矩，还是以BBH模型为例。
{:.info}
<!--more-->

# 模型及公式
BBH模型为

$$
\begin{aligned}
h^{q}(\mathbf{k})=& {\left[\gamma_{x}+\lambda_{x} \cos \left(k_{x}\right)\right] \Gamma_{4}+\lambda_{x} \sin \left(k_{x}\right) \Gamma_{3} } \\
&+\left[\gamma_{y}+\lambda_{y} \cos \left(k_{y}\right)\right] \Gamma_{2}+\lambda_{y} \sin \left(k_{y}\right) \Gamma_{1}
\end{aligned}
$$

实空间电四极矩计算公式为

$$
\begin{equation}
\begin{aligned}
q_{xy}&=\frac{1}{2\pi}\text{Im}\log[\text{det}(U^\dagger\hat{Q}U)\sqrt{\text{det}(Q^\dagger)}]\\
\hat{Q}\equiv\exp(i2\pi\hat{q}_{xy})\\
\hat{q}_{xy}\equiv\hat{x}\hat{y}/(L_xL_y)
\end{aligned}
\end{equation}
$$

这里$\hat{x}$和$\hat{y}$分别是沿着$x$和$y$方向的位置算符，而$L_x,L_y$则表示沿着此方向的长度。$U$是由占据态波函数构成的投影算符。且这里有下面的关系

$$
\sqrt{\text{det}(Q^\dagger)}=\exp(-i\pi\text{Tr}(\hat{q}_{xy}))
$$

那么利用上面的公式即可以在实空间中计算电四极矩。

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
    lamx::Float64 = 1.0  
    gamy::Float64 = gamx
    lamy::Float64 = lamx
    d0::Float64 = 0.00
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
#----------------------------------------
@everywhere function delta(x::Float64)
    gamma::Float64 = 0.01
    return 1.0/pi*gamma/(x*x + gamma*gamma)
end
#----------------------------------------
@everywhere function ldos()
    xn::Int64 = 30
    yn::Int64 = xn
    hn::Int64 = 4
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*hn
    omg::Float64 = 0.0
    val,vec = matset(xn,yn)
    # fx1 = "juldos-" * string(h0) * ".dat"
    fx1 = "dos.dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for ie in 1:N
                for i1 in 0:hn - 1
                    re1 = re1 + delta(val[ie] - omg)*(abs(vec[i0 + len2*i1,ie])^2)
                end
            end
            for ie in 1:Int(N/2)
                for i1 in 0:hn - 1
                    re2 += abs(vec[i0 + len2*i1,ie])^2 # 占据态波函数
                end
            end
            @printf(f1,"%5.2f\t%5.2f\t%15.8f\t%15.8f\n",ix,iy,re1,re2)
            # writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
end
#--------------------------------------------
@everywhere function Quadrupole()
    # 计算实空间电四极矩
    xn::Int64 = 30
    hn::Int64 = 4
    N::Int64 = xn^2*hn
    Nocc::Int64 = Int(N/2) # 占据态数量
    s0 = Diagonal(ones(hn,hn)) # 单位矩阵
    val,vec = matset(xn,xn,0.5)
    U = vec[:,1:Nocc] # 所有占据态波函数
    lx = Diagonal(range(1,xn,length = xn))
    ly = Diagonal(range(1,xn,length = xn))
    qxy = kron(kron(lx,ly)/(xn*xn),s0)
    Q = exp(2*pi*im*qxy)
    Qxy = angle(det(U'*Q*U))/(2.0*pi) - tr(qxy)/2.0
    return mod(Qxy,1.0)
end
#-----------------------------------------------
@everywhere function Quadrupole2(h0::Float64)
    # 计算实空间电四极矩
    xn::Int64 = 30
    hn::Int64 = 8
    N::Int64 = xn^2*hn
    Nocc::Int64 = Int(N/2) # 占据态数量
    s0 = Diagonal(ones(hn,hn)) # 单位矩阵
    val,vec = matset(xn,xn,h0)
    U = vec[:,1:Nocc] # 所有占据态波函数
    lx = Diagonal(range(1,xn,length = xn))
    ly = Diagonal(range(1,xn,length = xn))
    qxy = kron(kron(lx,ly)/(xn*xn),s0)
    Q = exp(2*pi*im*qxy)
    Qxy = angle(det(U'*Q*U))/(2.0*pi) - tr(qxy)/2.0
    return mod(Qxy,1.0)
end
#---------------------------------------------------
@everywhere function main2()
    h0list = -2:0.1:2
    re1 = SharedArray(zeros(Float64,length(h0list)))
    @sync @distributed for i0 in 1:length(h0list)
        re1[i0] = Quadrupole2(h0list[i0])
    end
    re1 = (a->(@sprintf "%15.8f" a)).(re1)
    h0list = (a->(@sprintf "%15.8f" a)).(h0list)
    f1 = open("Quadrupole-bbh.dat","w")
    writedlm(f1,[h0list re1])
    close(f1)
end
#---------------------------------------------------
@time main2()
```

这里顺手将代码写成了并行的，所以执行的时候为
```shell
julia -p NP code.jl
```
这里的`NP`就是描述想要用多少个核来记性计算，这里主要的计算过程并没有并行化，所以就只需要用`NP=1`即可。

# 参考

- 1.[Topological Phase Transitions in Disordered Electric Quadrupole Insulators](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.166801)


