---
title: 整理个自用的代码
tags: Julia 
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
```julia
using DelimitedFiles
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles
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
    writedlm(f1,map(real,val),"\t")
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
# main()
# charge(0.1)
@time wave(1.0)
```