---
title: 整理个自用的代码
tags: Julia Fortran
layout: article
license: true
toc: true
key: a20220521
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
最近稍微有点时间，把自己之前经常使用的一个代码做了一下精简，这里整理一下方便自己平时查阅使用。
{:.info}
<!--more-->
# 直接对角化
这里直接给出一个对角化的程序
```julia
using DelimitedFiles
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles,Printf
# --------------------------------------
@everywhere function boundary(xn::Int64,yn::Int64)
    len2::Int64 = xn*yn
    bry = zeros(Int64,8,len2)
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
            #-------------------------------------
            bry[5,i] = bry[1,i] + xn # right-above
            if iy == yn 
                bry[5,i] = bry[5,i] -len2
            end
            bry[6,i] = bry[2,i] + xn # left-above
            if iy == yn 
                bry[6,i] = bry[6,i] -len2
            end
            bry[7,i] = bry[1,i] - xn # right-below
            if iy == 1 
                bry[7,i] = bry[7,i] + len2
            end
            bry[8,i] = bry[2,i] - xn # left-below
            if iy == 1 
                bry[8,i] = bry[8,i] + len2
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
    g1 = kron(sz,s0,sz) # mass term
    g2 = kron(s0,sz,sx) # lambdax
    g3 = kron(sz,s0,sy) # lambday
    g4 = kron(sy,sy,s0) # dx^2-y^2
    g5 = kron(sx,sy,s0) # dxy
    g6 = kron(sz,sx,s0) # Zeeman
    g7 = kron(sz,s0,s0) # mu
    return g1,g2,g3,g4,g5,g6,g7
end
#---------------------------------------
@everywhere function matset(xn::Int64,yn::Int64,h0::Float64)
    m0::Float64 = 1.0
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    d0::Float64 = 0.0
    dx::Float64 = 0.5
    dy::Float64 = -dx
    dp::Float64 = 0.3
    mu::Float64 = 0.1
    # h0::Float64 = 1.0
    hn::Int64 = 8
    N::Int64 = xn*yn*hn
    len2::Int64 = xn*yn
    ham = zeros(ComplexF64,N,N)
    g1,g2,g3,g4,g5,g6,g7 = gamma()
    bry = boundary(xn,yn)
    #---
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            for i1 in 0:hn -1 
                for i2 in 0:hn - 1
                    ham[i0 + len2*i1,i0 + len2*i2] = m0*g1[i1 + 1,i2 + 1] + d0*g4[i1 + 1,i2 + 1] - mu*g7[i1 + 1,i2 + 1] + h0*g6[i1 + 1,i2 + 1]
                    if ix != xn
                        ham[i0 + len2*i1,bry[1,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] + ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if ix != 1
                        ham[i0 + len2*i1,bry[2,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] - ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if iy != yn
                        ham[i0 + len2*i1,bry[3,i0] + len2*i2] = -ty/2.0*g1[i1 + 1,i2 + 1] + ay/(2*im)*g3[i1 + 1,i2 + 1] + dy/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if iy != 1
                        ham[i0 + len2*i1,bry[4,i0] + len2*i2] = -ty/2.0*g1[i1 + 1,i2 + 1] - ay/(2*im)*g3[i1 + 1,i2 + 1] + dy/2.0*g4[i1 + 1,i2 + 1]
                    end
                end
            end
        end
    end
    #-----------------
    # for m1 in 1:N
    #     for m2 in 1:N
    #         if ham[m1,m2] != conj(ham[m2,m1])
    #             println("Hamiltonian is not hermitian")
    #             # println("(",m1,",",m2,")",ham[m1,m2],ham[m2,m1])
    #             break
    #         end
    #     end
    # end
    #-----------------------------------------
    if ishermitian(ham)
        val,vec = eigen(ham)
    else
        println("Hamiltonian is not hermitian")
        # break
    end
    fx1 = "juliaval-" * string(h0) * ".dat"
    f1 = open(fx1,"w")
    # writedlm(f1,map(real,val),"\t")
    ind = (a->(@sprintf "%15.8f" a)).(range(1,length(val),length = length(val))) # 本征值格式化输出
    val2 = (a->(@sprintf "%15.8f" a)).(map(real,val))
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
@everywhere function ldos(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "juldos-" * string(h0) * ".dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for ie in 1:N
                for i1 in 0:7
                    re1 = re1 + delta(val[ie] - omg)*(abs(vec[i0 + len2*i1,ie])^2)
                end
            end
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            @printf(f1,"%15.8f\t%15.8f\t%15.8f\t%15.8f\n",ix,iy,re1,re2)
            # writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
end
#-----------------------------------------
@everywhere function wave(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "wave-" * string(h0) * "dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for m1 = 0:7
                for m2 = -1:2
                    re1 = re1 + abs(vec[i0 + m1*len2, Int(N/2) + m2])^2
                end 
            end 
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            @printf(f1,"%15.8f\t%15.8f\t%15.8f\t%15.8f\n",ix,iy,re1,re2)
            # writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
    
end
#----------------------------------------
@everywhere function main()
    @sync @distributed for h0 in -1:0.05:1
        ldos(h0)
        # charge(h0)
        # println(h0)
    end
end
#--------------------------------------
# main()
# charge(0.1)
@time wave(1.0)
```

# 限定范围本征值
通常在问题研究的时候，其实并不需要所有的本征值，这里通过`eigen`的参数设置只求解确定范围内的本征值
```julia
using DelimitedFiles
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles,Printf
# --------------------------------------
@everywhere function boundary(xn::Int64,yn::Int64)
    len2::Int64 = xn*yn
    bry = zeros(Int64,8,len2)
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
            #-------------------------------------
            bry[5,i] = bry[1,i] + xn # right-above
            if iy == yn 
                bry[5,i] = bry[5,i] -len2
            end
            bry[6,i] = bry[2,i] + xn # left-above
            if iy == yn 
                bry[6,i] = bry[6,i] -len2
            end
            bry[7,i] = bry[1,i] - xn # right-below
            if iy == 1 
                bry[7,i] = bry[7,i] + len2
            end
            bry[8,i] = bry[2,i] - xn # left-below
            if iy == 1 
                bry[8,i] = bry[8,i] + len2
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
    g1 = kron(sz,s0,sz) # mass term
    g2 = kron(s0,sz,sx) # lambdax
    g3 = kron(sz,s0,sy) # lambday
    g4 = kron(sy,sy,s0) # dx^2-y^2
    g5 = kron(sx,sy,s0) # dxy
    g6 = kron(sz,sx,s0) # Zeeman
    g7 = kron(sz,s0,s0) # mu
    return g1,g2,g3,g4,g5,g6,g7
end
#---------------------------------------
@everywhere function matset(xn::Int64,yn::Int64,h0::Float64)
    m0::Float64 = 1.0
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    d0::Float64 = 0.0
    dx::Float64 = 0.5
    dy::Float64 = -dx
    dp::Float64 = 0.3
    mu::Float64 = 0.1
    # h0::Float64 = 1.0
    hn::Int64 = 8
    N::Int64 = xn*yn*hn
    len2::Int64 = xn*yn
    ham = zeros(ComplexF64,N,N)
    g1,g2,g3,g4,g5,g6,g7 = gamma()
    bry = boundary(xn,yn)
    #---
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            for i1 in 0:hn -1 
                for i2 in 0:hn - 1
                    ham[i0 + len2*i1,i0 + len2*i2] = m0*g1[i1 + 1,i2 + 1] + d0*g4[i1 + 1,i2 + 1] - mu*g7[i1 + 1,i2 + 1] + h0*g6[i1 + 1,i2 + 1]
                    if ix != xn
                        ham[i0 + len2*i1,bry[1,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] + ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if ix != 1
                        ham[i0 + len2*i1,bry[2,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] - ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if iy != yn
                        ham[i0 + len2*i1,bry[3,i0] + len2*i2] = -ty/2.0*g1[i1 + 1,i2 + 1] + ay/(2*im)*g3[i1 + 1,i2 + 1] + dy/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if iy != 1
                        ham[i0 + len2*i1,bry[4,i0] + len2*i2] = -ty/2.0*g1[i1 + 1,i2 + 1] - ay/(2*im)*g3[i1 + 1,i2 + 1] + dy/2.0*g4[i1 + 1,i2 + 1]
                    end
                end
            end
        end
    end
    #-----------------
    # for m1 in 1:N
    #     for m2 in 1:N
    #         if ham[m1,m2] != conj(ham[m2,m1])
    #             println("Hamiltonian is not hermitian")
    #             # println("(",m1,",",m2,")",ham[m1,m2],ham[m2,m1])
    #             break
    #         end
    #     end
    # end
    #-----------------------------------------
    # if ishermitian(ham)
    #     val,vec = eigen(ham)
    # else
    #     println("Hamiltonian is not hermitian")
    #     # break
    # end
    # fx1 = "juliaval-" * string(h0) * ".dat"
    # f1 = open(fx1,"w")
    # writedlm(f1,map(real,val),"\t")
    # close(f1)
    # return map(real,val),vec
    return Hermitian(ham)
end
#----------------------------------------
@everywhere function delta(x::Float64)
    gamma::Float64 = 0.01
    return 1.0/pi*gamma/(x*x + gamma*gamma)
end
#----------------------------------------
@everywhere function ldos(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "juldos-" * string(h0) * ".dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for ie in 1:N
                for i1 in 0:7
                    re1 = re1 + delta(val[ie] - omg)*(abs(vec[i0 + len2*i1,ie])^2)
                end
            end
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            @printf(f1,"%15.8f\t%15.8f\t%15.8f\t%15.8f\n",ix,iy,re1,re2)
            # writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
end
#-----------------------------------------
@everywhere function wave(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "wave-" * string(h0) * "dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for m1 = 0:7
                for m2 = -1:2
                    re1 = re1 + abs(vec[i0 + m1*len2, Int(N/2) + m2])^2
                end 
            end 
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            @printf(f1,"%15.8f\t%15.8f\t%15.8f\t%15.8f\n",ix,iy,re1,re2)
            # writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
    
end
#----------------------------------------
@everywhere function main()
    @sync @distributed for h0 in -1:0.05:1
        ldos(h0)
        # charge(h0)
        # println(h0)
    end
end
#--------------------------------------
@everywhere function test(xn::Int64)
    ham = matset(xn,xn,1.0)
    val = eigvals(ham,-1.0,1.0)
    #val = eigvals(ham)
    f1 = open("val-1.dat","w")
    writedlm(f1,val)
    close(f1)
end 
#--------------------------------------
@sync @distributed for i0 in 10:50
    @time test(Int(i0))
end
```

这里主要是调整了函数参数
```julia
val = eigvals(ham,-1.0,1.0) 
```
将本征值求解限定在了`[-1,1]`这个区间内，而且在构造矩阵的时候必须要使用`Hermitian`这个函数让矩阵是厄米的才可以这样做。通过测试发现这种限制确定范围内的本征值的方法并不能在求解速度上带来提升，所以只能是方便自己在处理特定问题的时候，只是输出确定范围内的本征值。

# 稀疏矩阵对角化
上面用到的方法虽然可以给出制定范围内的本征值，但是求解速度可能并不是最优的，这里进一步用`Arpack`这个包中的函数`eigs`可以更加方便快速的得到一定数量的本征值和本征矢量
```julia
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles,Arpack
# --------------------------------------
@everywhere function boundary(xn::Int64,yn::Int64)
    len2::Int64 = xn*yn
    bry = zeros(Int64,8,len2)
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
            #-------------------------------------
            bry[5,i] = bry[1,i] + xn # right-above
            if iy == yn 
                bry[5,i] = bry[5,i] -len2
            end
            bry[6,i] = bry[2,i] + xn # left-above
            if iy == yn 
                bry[6,i] = bry[6,i] -len2
            end
            bry[7,i] = bry[1,i] - xn # right-below
            if iy == 1 
                bry[7,i] = bry[7,i] + len2
            end
            bry[8,i] = bry[2,i] - xn # left-below
            if iy == 1 
                bry[8,i] = bry[8,i] + len2
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
    g1 = kron(sz,s0,sz) # mass term
    g2 = kron(s0,sz,sx) # lambdax
    g3 = kron(sz,s0,sy) # lambday
    g4 = kron(sy,sy,s0) # dx^2-y^2
    g5 = kron(sx,sy,s0) # dxy
    g6 = kron(sz,sx,s0) # Zeeman
    g7 = kron(sz,s0,s0) # mu
    return g1,g2,g3,g4,g5,g6,g7
end
#---------------------------------------
@everywhere function matset(xn::Int64,yn::Int64,h0::Float64)
    m0::Float64 = 1.0
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    d0::Float64 = 0.0
    dx::Float64 = 0.5
    dy::Float64 = -dx
    dp::Float64 = 0.3
    mu::Float64 = 0.1
    # h0::Float64 = 1.0
    hn::Int64 = 8
    N::Int64 = xn*yn*hn
    len2::Int64 = xn*yn
    ham = zeros(ComplexF64,N,N)
    g1,g2,g3,g4,g5,g6,g7 = gamma()
    bry = boundary(xn,yn)
    #---
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            for i1 in 0:hn -1 
                for i2 in 0:hn - 1
                    ham[i0 + len2*i1,i0 + len2*i2] = m0*g1[i1 + 1,i2 + 1] + d0*g4[i1 + 1,i2 + 1] - mu*g7[i1 + 1,i2 + 1] + h0*g6[i1 + 1,i2 + 1]
                    if ix != xn
                        ham[i0 + len2*i1,bry[1,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] + ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if ix != 1
                        ham[i0 + len2*i1,bry[2,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] - ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if iy != yn
                        ham[i0 + len2*i1,bry[3,i0] + len2*i2] = -ty/2.0*g1[i1 + 1,i2 + 1] + ay/(2*im)*g3[i1 + 1,i2 + 1] + dy/2.0*g4[i1 + 1,i2 + 1]
                    end
                end
            end
        end
    end
    #-----------------
    # for m1 in 1:N
    #     for m2 in 1:N
    #         if ham[m1,m2] != conj(ham[m2,m1])
    #             println("Hamiltonian is not hermitian")
    #             # println("(",m1,",",m2,")",ham[m1,m2],ham[m2,m1])
    #             break
    #         end
    #     end
    # end
    #-----------------------------------------
    # if ishermitian(ham)
    #     val,vec = eigen(ham)
    # else
    #     println("Hamiltonian is not hermitian")
    #     # break
    # end
    # fx1 = "juliaval-" * string(h0) * ".dat"
    # f1 = open(fx1,"w")
    # writedlm(f1,map(real,val),"\t")
    # close(f1)
    # return map(real,val),vec
    return Hermitian(ham)
end
#----------------------------------------
@everywhere function delta(x::Float64)
    gamma::Float64 = 0.01
    return 1.0/pi*gamma/(x*x + gamma*gamma)
end
#----------------------------------------
@everywhere function ldos(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "juldos-" * string(h0) * ".dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for ie in 1:N
                for i1 in 0:7
                    re1 = re1 + delta(val[ie] - omg)*(abs(vec[i0 + len2*i1,ie])^2)
                end
            end
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
end
#-----------------------------------------
@everywhere function wave(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "wave-" * string(h0) * "dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for m1 = 0:7
                for m2 = -1:2
                    re1 = re1 + abs(vec[i0 + m1*len2, Int(N/2) + m2])^2
                end 
            end 
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
    
end
#----------------------------------------
@everywhere function main()
    @sync @distributed for h0 in -1:0.05:1
        ldos(h0)
        # charge(h0)
        # println(h0)
    end
end
#--------------------------------------
@everywhere function test(xn::Int64)
    ham = matset(xn,xn,1.0)
    val,vec = eigs(ham,nev = 100,maxiter=1500,which=:SM) 
    f1 = open("val-3.dat","w")
    writedlm(f1,sort(map(real,val)))
    close(f1)
end 
#--------------------------------------
f1 = open("time-2.dat","w")
for i0 in 1:50
    t = @elapsed test(30)
    writedlm(f1,t)
end
close(f1)
```
这里最主要的就是利用
```julia
val,vec = eigs(ham,nev = 100,maxiter=1500,which=:SM) 
```
这里`nev`设置想要多少个本征值，`maxiter`设置算法迭代的次数，`which`可以设置想要求解最小`SM`的`nev`个本征值还是最大`LM`的`nev`个本征值。

这个方法只在获取少量本征值的时候速度较快，如果获取本征值的数量就多，就直接使用`eigen`来操作，我测试了一下，同样是30*30的格点数目，7200维的矩阵对角化，通过`eigs`获取所有的本征值，耗时非常长，远大于`eigen`所需要的时长。
{:.warning}

# 格式化输出

前面计算的时候，输出的数据格式都是不整齐的，这里将代码优化一下，让输出的数据更加整齐一下
```julia
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles,Arpack,Printf
# --------------------------------------
@everywhere function boundary(xn::Int64,yn::Int64)
    len2::Int64 = xn*yn
    bry = zeros(Int64,8,len2)
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
            #-------------------------------------
            bry[5,i] = bry[1,i] + xn # right-above
            if iy == yn 
                bry[5,i] = bry[5,i] -len2
            end
            bry[6,i] = bry[2,i] + xn # left-above
            if iy == yn 
                bry[6,i] = bry[6,i] -len2
            end
            bry[7,i] = bry[1,i] - xn # right-below
            if iy == 1 
                bry[7,i] = bry[7,i] + len2
            end
            bry[8,i] = bry[2,i] - xn # left-below
            if iy == 1 
                bry[8,i] = bry[8,i] + len2
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
    g1 = kron(sz,s0,sz) # mass term
    g2 = kron(s0,sz,sx) # lambdax
    g3 = kron(sz,s0,sy) # lambday
    g4 = kron(sy,sy,s0) # dx^2-y^2
    g5 = kron(sx,sy,s0) # dxy
    g6 = kron(sz,sx,s0) # Zeeman
    g7 = kron(sz,s0,s0) # mu
    return g1,g2,g3,g4,g5,g6,g7
end
#---------------------------------------
@everywhere function matset(xn::Int64,yn::Int64,h0::Float64)
    m0::Float64 = 1.0
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    d0::Float64 = 0.0
    dx::Float64 = 0.5
    dy::Float64 = -dx
    dp::Float64 = 0.3
    mu::Float64 = 0.1
    # h0::Float64 = 1.0
    hn::Int64 = 8
    N::Int64 = xn*yn*hn
    len2::Int64 = xn*yn
    ham = zeros(ComplexF64,N,N)
    g1,g2,g3,g4,g5,g6,g7 = gamma()
    bry = boundary(xn,yn)
    #---
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            for i1 in 0:hn -1 
                for i2 in 0:hn - 1
                    ham[i0 + len2*i1,i0 + len2*i2] = m0*g1[i1 + 1,i2 + 1] + d0*g4[i1 + 1,i2 + 1] - mu*g7[i1 + 1,i2 + 1] + h0*g6[i1 + 1,i2 + 1]
                    if ix != xn
                        ham[i0 + len2*i1,bry[1,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] + ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if ix != 1
                        ham[i0 + len2*i1,bry[2,i0] + len2*i2] = -tx/2.0*g1[i1 + 1,i2 + 1] - ax/(2*im)*g2[i1 + 1,i2 + 1] + dx/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if iy != yn
                        ham[i0 + len2*i1,bry[3,i0] + len2*i2] = -ty/2.0*g1[i1 + 1,i2 + 1] + ay/(2*im)*g3[i1 + 1,i2 + 1] + dy/2.0*g4[i1 + 1,i2 + 1]
                    end
                    if iy != 1
                        ham[i0 + len2*i1,bry[4,i0] + len2*i2] = -ty/2.0*g1[i1 + 1,i2 + 1] - ay/(2*im)*g3[i1 + 1,i2 + 1] + dy/2.0*g4[i1 + 1,i2 + 1]
                    end
                end
            end
        end
    end
    #-----------------
    # for m1 in 1:N
    #     for m2 in 1:N
    #         if ham[m1,m2] != conj(ham[m2,m1])
    #             println("Hamiltonian is not hermitian")
    #             # println("(",m1,",",m2,")",ham[m1,m2],ham[m2,m1])
    #             break
    #         end
    #     end
    # end
    #-----------------------------------------
    # if ishermitian(ham)
    #     val,vec = eigen(ham)
    # else
    #     println("Hamiltonian is not hermitian")
    #     # break
    # end
    # fx1 = "juliaval-" * string(h0) * ".dat"
    # f1 = open(fx1,"w")
    # writedlm(f1,map(real,val),"\t")
    # close(f1)
    # return map(real,val),vec
    return Hermitian(ham)
end
#----------------------------------------
@everywhere function delta(x::Float64)
    gamma::Float64 = 0.01
    return 1.0/pi*gamma/(x*x + gamma*gamma)
end
#----------------------------------------
@everywhere function ldos(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "juldos-" * string(h0) * ".dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for ie in 1:N
                for i1 in 0:7
                    re1 = re1 + delta(val[ie] - omg)*(abs(vec[i0 + len2*i1,ie])^2)
                end
            end
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
end
#-----------------------------------------
@everywhere function wave(h0::Float64)
    xn::Int64 = 10
    yn::Int64 = xn
    len2::Int64 = xn*yn
    N::Int64 = xn*yn*8
    omg::Float64 = 0.0
    val,vec = matset(xn,yn,h0)
    fx1 = "wave-" * string(h0) * "dat"
    f1 = open(fx1,"w")
    for iy in 1:yn
        for ix in 1:xn
            i0 = (iy - 1)*xn + ix
            re1 = 0
            re2 = 0
            for m1 = 0:7
                for m2 = -1:2
                    re1 = re1 + abs(vec[i0 + m1*len2, Int(N/2) + m2])^2
                end 
            end 
            for ie in 1:Int(N/2)
                for i1 in 0:7
                    re2 += abs(vec[i0 + len2*i1,ie])^2
                end
            end
            writedlm(f1,[ix iy re1 re2],"\t")
        end
    end
    close(f1)
end
#----------------------------------------
@everywhere function main()
    @sync @distributed for h0 in -1:0.05:1
        ldos(h0)
        # charge(h0)
        # println(h0)
    end
end
#--------------------------------------
@everywhere function test(xn::Int64,valnum::Int64)
    ham = matset(xn,xn,1.0)
    val,vec = eigs(ham,nev = valnum,maxiter=1500,which=:SM) # 该函数可以获得确定范围内的本征值，而且计算速度提升较大
    f1 = open("val.dat","w")
    f2 = open("val-order.dat","w")
    writedlm(f1,map(real,val))
    writedlm(f2,sort(map(real,val)))
    close(f1)
end 
#----------------------------------------
@everywhere function test()
    # 测试对角化所需要的时间，这里固定了格点大小
    f1 = open("time-size.dat","w")
    for i0 in 1:50
        t = @elapsed test(30,100) #格点大小为30*30，求解最小的100个本征值
        writedlm(f1,t)
    end
    close(f1)
end
#--------------------------------------
@everywhere function test2()
    # 测试对角化的到不同数量的本征值所需要的时间，这里固定了格点大小
    f1 = open("time-valnum.dat","w")
    for i0 in 100:100:7200
        t = @elapsed test(30,Int(i0)) #格点大小为30*30，求解最小的100个本征值
        writedlm(f1,[i0 t])
    end
    close(f1)
end
#---------------------------------------
@everywhere function main1()
    xn::Int64 = 10
    valnum::Int64 = 100
    ham = matset(xn,xn,1.0)
    val,vec = eigs(ham,nev = valnum,maxiter=1500,which=:SM) # 该函数可以获得确定范围内的本征值，而且计算速度提升较大    
    f1 = open("val-format.dat","w")
    val =  map(real,val)
    te1 = sortperm(val) # 对本征值给个索引顺序
    val =  (a->(@sprintf "%15.8g" a)).(map(real,val)) 
    te2 =  (a->(@sprintf "%15.8g" a)).(te1)
    writedlm(f1,[val  te2 val[te1]])
    close(f1)
end
#--------------------------------------
@time main1()
```
这里主要修改的就是
```julia
val =  (a->(@sprintf "%15.8g" a)).(map(real,val)) 
te2 =  (a->(@sprintf "%15.8g" a)).(te1)
```
通过打印的方式将数据整理成整齐的格式，这样方便查看和后续对数据的处理。同样还可以将所有的数据都处理成浮点类型，保持整洁性
```julia
val =  (a->(@sprintf "%15.8f" a)).(map(real,val)) 
te2 =  (a->(@sprintf "%15.8f" a)).(te1)
```
除了上面的方式之外，还可以使用这些函数直接将输出定向到文件中。
```julia
function test3()
    xn::Int64 = 1000
    f1 = open("format.dat","w")
    for i0 in 1:xn
        @printf(f1,"%15.8f\t %15.8f\n",i0,sin(i0))
    end
    close(f1)
end
```
这样即可以实现文件的格式化输出的，自己暂时也算是解决的在`julia`环境下数据的格式化操作。

# Fortran 版本
在有些计算上还是`Fortran`更快一些，所以这里也把这个方法用`Fotran`重写了一下，因为自己之前的好多计算也都是用`Fortran`写的
```fortran
    ! Author:YuXuanLi
    ! E-Mail:yxli406@gmail.com
    ! 随空间位置变化的微分电导谱
    module pub
    implicit none
    integer xn,yn,N,len2,hn,omgN!(计算ldos时撒点数量)
    parameter(xn = 30,yn = xn,hn = 8,N = xn*yn*hn,len2 = xn*yn, omgN = 100)
    complex,parameter::im = (0.0,1.0) 
    real,parameter::pi = 3.14159265359
    complex Ham(N,N) 
    integer bry(8,len2)
    real m0,omega,mu
    real tx,ty,del 
    real d0,dx,dy,dp  
    real ax,ay,h0
    complex g1(hn,hn),g2(hn,hn),g3(hn,hn),g4(hn,hn)
    complex g5(hn,hn),g6(hn,hn),g7(hn,hn),g8(hn,hn)
    !-------------------lapack parameter----------------------------------------
    integer::lda = N
    integer,parameter::lwmax=2*N+N**2
    real,allocatable::w(:)
    complex*8,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   
    integer lrwork    
    integer liwork   
    integer info
    end module pub
!========== PROGRAM START ==========================
    program sol
    use pub
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1 + 5*N + 2*N**2))
    allocate(iwork(3 + 5*N))
    !-------------------
    m0 = 1.
    tx = 2.0   
    ty = 2.0  
    ! mu = -0.3 
    ax = 2.0
    ay = 2.0 
    d0 = 0.0
    dx = 0.5
    dy = -dx
    mu = 0.
    h0 = 1.0 ! 这个参数下面有corner态
    dp = 0.3
    call matset(1)
    call ldos(1)
    stop
    end program sol
!===========================Local Density of State=============================
    subroutine ldos(m2)
    ! 计算确定能量omg下面所有位置的ldos
    use pub
    integer m,l1,k,l2,m2
    real omg,re1,re2,re3
    character*20::str1,str2,str3,str4
    real,external::delta
    omg = 0.0
    str1 = "ldos-"
    write(str2,"(I4.4)")m2
    str3 = ".dat"
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(12,file=str4) 
    do l1 = 1,yn
        do l2 = 1,xn
            k = l2 + (l1-1)*xn
            re1 = 0
            do m=1,N
                do i1 = 0,hn - 1
                   re1 = re1 + delta(w(m) - omg)*abs(Ham(k + len2 * i1, m))**2 ! ldos
                end do
            end do
            re2 = 0
            do m = N/2 - 2,N/2 + 1
                do i1 = 0,hn - 1
                    re2 = re2 + abs(Ham(k + len2 * i1, m))**2 ! 零能波函数分布
                end do
            end do
            re3 = 0
            do m = 1,N/2
                do i1 = 0,hn - 1
                    re3 = re3 + abs(Ham(k + len2 * i1, m))**2 ! 占据态求和
                end do
            end do
            write(12,999)1.*l1,1.*l2,re1,re2,re3
        end do
    end do
    close(12)
    999 format(10f11.6)
    return
    end subroutine ldos
!======================================================================
    real function delta(x)
    real x
    real::gamma = 0.01
    delta = 1.0/3.1415926535*gamma/(x*x + gamma*gamma)
    end function delta
!=======================================================
    subroutine matset(input)
    use pub
    integer input
    integer i1,i2,ix,iy,i0
    call pauli()
    call boundary()
    do iy = 1,yn
        do ix = 1,xn
            i0 = (iy - 1)*xn + ix
            do i1 = 0,hn -1 
                do i2 = 0,hn - 1
                    ham(i0 + len2 * i1,i0 + len2 * i2) = m0*g1(i1 + 1,i2 + 1) + d0*g4(i1 + 1,i2 + 1) - mu*g7(i1 + 1,i2 + 1) + h0*g5(i1 + 1,i2 + 1)
                    if(ix.ne.xn)ham(i0 + len2 * i1,bry(1,i0) + len2 * i2) = -tx/2.0*g1(i1 + 1,i2 + 1) + ax/(2*im)*g2(i1 + 1,i2 + 1) + dx/2.0*g4(i1 + 1,i2 + 1)
                    if(ix.ne.1)ham(i0 + len2 * i1,bry(2,i0) + len2 * i2) = -tx/2.0*g1(i1 + 1,i2 + 1) - ax/(2*im)*g2(i1 + 1,i2 + 1) + dx/2.0*g4(i1 + 1,i2 + 1)
                    if(iy.ne.yn)ham(i0 + len2 * i1,bry(3,i0) + len2 * i2) = -ty/2.0*g1(i1 + 1,i2 + 1) + ay/(2*im)*g3(i1 + 1,i2 + 1) + dy/2.0*g4(i1 + 1,i2 + 1)
                    if(iy.ne.1)ham(i0 + len2 * i1,bry(4,i0) + len2 * i2) = -ty/2.0*g1(i1 + 1,i2 + 1) - ay/(2*im)*g3(i1 + 1,i2 + 1) + dy/2.0*g4(i1 + 1,i2 + 1)
                end do
            end do
        end do
    end do
    call isHermitian()
    call eigsol(input)
    return
    end subroutine matset
!==========================================================
    subroutine pauli()
    use pub
    !---------Kinetic energy
    g1(1,1) = 1
    g1(2,2) = -1
    g1(3,3) = 1
    g1(4,4) = -1
    g1(5,5) = -1
    g1(6,6) = 1
    g1(7,7) = -1
    g1(8,8) = 1
    !----------SOC-x
    g2(1,2) = 1
    g2(2,1) = 1
    g2(3,4) = -1
    g2(4,3) = -1
    g2(5,6) = 1
    g2(6,5) = 1
    g2(7,8) = -1
    g2(8,7) = -1
    !---------SOC-y
    g3(1,2) = -im
    g3(2,1) = im
    g3(3,4) = -im
    g3(4,3) = im
    g3(5,6) = im 
    g3(6,5) = -im
    g3(7,8) = im
    g3(8,7) = -im
    !------------------- dx^2-y^2
    g4(1,7) = -1
    g4(2,8) = -1
    g4(3,5) = 1
    g4(4,6) = 1
    g4(7,1) = -1
    g4(8,2) = -1
    g4(5,3) = 1
    g4(6,4) = 1
    !-------------------- layer couple
    g5(1,3) = 1
    g5(2,4) = 1
    g5(3,1) = 1
    g5(4,2) = 1
    g5(5,7) = -1
    g5(6,8) = -1
    g5(7,5) = -1
    g5(8,6) = -1
    !-------------------- dxy pairing
    g6(1,7) = -im
    g6(2,8) = -im
    g6(3,5) = im
    g6(4,6) = im
    g6(5,3) = -im
    g6(6,4) = -im
    g6(7,1) = im
    g6(8,2) = im
    !-------------------- dxy pairing
    g7(1,1) = 1
    g7(2,2) = 1
    g7(3,3) = 1
    g7(4,4) = 1
    g7(5,5) = -1
    g7(6,6) = -1
    g7(7,7) = -1
    g7(8,8) = -1
    return
    end subroutine pauli
!===========================================================
    subroutine boundary()
    use pub
    integer i,ix,iy
    bry = 0
    do iy = 1,yn  
        do ix = 1,xn 
            i = (iy-1)*xn + ix
            bry(1,i) = i + 1    
            if(ix.eq.xn)bry(1,i) = bry(1,i) - xn    
            bry(2,i) = i - 1    
            if(ix.eq.1)bry(2,i) = bry(2,i) + xn     
            bry(3,i) = i + xn   
            if(iy.eq.yn)bry(3,i) = bry(3,i) - len2  
            bry(4,i)= i - xn    
            if(iy.eq.1)bry(4,i) = bry(4,i) + len2   
            !--------NNN-----------------------------
            bry(5,i) = bry(1,i) + xn    ! right-above
            if(iy.eq.yn)bry(5,i) = bry(5,i) -len2
            bry(6,i) = bry(2,i) + xn    ! left-above
            if(iy.eq.yn)bry(6,i) = bry(6,i) -len2
            bry(7,i) = bry(1,i) - xn    ! right-below
            if(iy.eq.1)bry(7,i) = bry(7,i) + len2
            bry(8,i) = bry(2,i) - xn    ! left-below
            if(iy.eq.1)bry(8,i) = bry(8,i) + len2
        enddo
    enddo
    return
    end subroutine boundary
!============================================================
    subroutine isHermitian()
    use pub
    integer i,j
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                open(16,file = 'hermitian.dat')
                write(16,*)i,j
                write(16,*)Ham(i,j)
                write(16,*)Ham(j,i)
                write(*,*)"Hamiltonian is not Hermitian"
                stop
            end if
        end do
    end do
    close(16)
    return
    end subroutine isHermitian
!================= Hermita= Matrices solve ==============
    subroutine eigsol(input)
    use pub
    integer m,input
    character*20::str1,str2,str3,str4
    str1 = "eigval-"
    str3 = ".dat"
    write(str2,"(I4.4)")input
    str4 = trim(str1)//trim(str2)//trim(str3)
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','U',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','U',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(110,file="mes.dat",status="unknown")
        write(110,*)'The algorithm failed to compute eigenvalues.'
        close(110)
    end if
    open(120,file = str4) 
    do m = 1,N
        write(120,998)1.0 * m,w(m)
    end do
    close(120)
998 format(10f11.6)
    return
    end subroutine eigsol
```