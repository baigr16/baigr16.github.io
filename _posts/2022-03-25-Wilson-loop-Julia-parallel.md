---
title: Julia 并行计算Wilson loop
tags: MPI Julia 
layout: article
license: true
toc: true
key: a202203025
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

这里使用julia并行计算Wilson loop。
{:.info}
<!--more-->

# 前言
这里还是使用BHZ模型为例来做计算，其实主要也就是稍微将代码进行了一些改动，使得其能进行并行计算，提高计算效率。
# 代码
```julia
using DelimitedFiles
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles
@everywhere function pauli()
    hn::Int64 = 2
    g0 = zeros(ComplexF64,hn,hn)
    gx = zeros(ComplexF64,hn,hn)
    gy = zeros(ComplexF64,hn,hn)
    gz = zeros(ComplexF64,hn,hn)
    #---------------
    g0[1,1] = 1
    g0[2,2] = 1

    gx[1,2] = 1
    gx[2,1] = 1

    gy[1,2] = -im
    gy[2,1] = im

    gz[1,1] = 1
    gz[2,2] = -1
    return g0,gx,gy,gz
end
#---------------------------------------------------------------
@everywhere function hamset(kx::Float64,ky::Float64,m0::Float64)::Matrix{ComplexF64}
    # m0::Float64 = 1.5
    tx::Float64 = 1.
    ty::Float64 = 1.
    ax::Float64 = 1.
    ay::Float64 = 1.
    s0,sx,sy,sz = pauli()
    ham = (m0 - tx*cos(kx) - ty*cos(ky))*kron(s0,sz) + ax*sin(kx)*kron(sz,sx) + ay*sin(ky)*kron(s0,sy) 
    return ham
end
#--------------------------------------------------------------------------------------
@everywhere function Wilson_kx(kx::Float64,m0::Float64)
    nky::Int64 = 100
    hn::Int64 = 4
    Nocc::Int64 = Int(hn/2)
    wave = zeros(ComplexF64,hn,hn,nky) # 存储哈密顿量对应的波函数
    wan = zeros(ComplexF64,Nocc,Nocc) # 存储Wannier哈密顿量对应的波函数
    F = zeros(ComplexF64,Nocc,Nocc)
    vec2 = zeros(ComplexF64,Nocc,Nocc)  # 重新排列顺序后Wannier Hamiltonian的本征矢量
    kylist = range(0, 2*pi, length = nky)
    for iy in 1:nky # 固定kx，沿着ky方向计算Wilson loop
        ky = kylist[iy]
        val,vec = eigen(hamset(kx,ky,m0))
        wave[:,:,iy] = vec[:,:] # 存储波函数
    end
    wave[:,:,nky] = wave[:,:,1] # 在边界上波函数首尾相接
    for i1 in 1:Nocc
        F[i1,i1] = 1 # 构建单位矩阵
    end
    for i1 in 1:nky - 1 # index ki lattice
        for i2 in 1:Nocc
            for i3 in 1:Nocc
                wan[i2,i3] = wave[:,i2,i1]' * wave[:,i3,i1 + 1]   # 计算Berry联络
            end
        end
        F = wan * F # 沿着ky方向构造Wannier哈密顿量
    end
    val,vec = eigen(F) # 对求解得到的本征矢量按照本征值大小排列
    val = map(angle,val)
    for i0 in 1:length(val)
        if val[i0] < 0
            val[i0] += 1
        end
    end
    return val
end
#-------------------------------------------------------------------------------------
@everywhere function Wilson_ky(ky::Float64,m0::Float64)
    nkx::Int64 = 100
    hn::Int64 = 4
    Nocc::Int64 = Int(hn/2)
    wave = zeros(ComplexF64,hn,hn,nkx)
    wan = zeros(ComplexF64,Nocc,Nocc)
    F = zeros(ComplexF64,Nocc,Nocc)
    vec2 = zeros(ComplexF64,Nocc,Nocc)  # 重新排列顺序后Wannier Hamiltonian的本征矢量
    kxlist = range(0, 2*pi, length = nkx)
    for ix in 1:nkx # 固定ky，沿着kx方向计算Wilson loop
        kx = kxlist[ix]
        val,vec = eigen(hamset(kx,ky,m0))
        wave[:,:,ix] = vec[:,:]
    end
    wave[:,:,nkx] = wave[:,:,1] # 波函数首尾相接
    for i1 in 1:Nocc
        F[i1,i1] = 1
    end
    for i1 in 1:nkx - 1
        for i2 in 1:Nocc
            for i3 in 1:Nocc
                wan[i2,i3] = wave[:,i2,i1]' * wave[:,i3,i1 + 1]   # 计算Berry联络
            end
        end
        F = wan * F # 构造Wannier 哈密顿量
    end
    val,vec = eigen(F) # 对求解得到的本征矢量按照本征值大小排列
    val = map(angle,val)
    for i0 in 1:length(val)
        if val[i0] < 0
            val[i0] += 1
        end
    end
    return val  
end
#-----------------------------------------------------------------------------------------------
@everywhere function Wilsonloop(m0::Float64,cont::Int64)
    # 改变质量看Wilson loop变化
    klist = -pi:0.1:pi
    re1 = zeros(Float64,length(klist),2)
    re2 = zeros(Float64,length(klist),2)
    i0 = 0
    for kx in klist
        i0 += 1
        x1 = Wilson_kx(kx,m0)
        x2 = Wilson_ky(kx,m0)
        re1[i0,:] = x1
        re2[i0,:] = x2
    end
    fx1 = "Wilson-kx-" * string(cont) * ".dat"
    f1 = open(fx1,"w")
    writedlm(f1,[klist re1],"\t")
    close(f1)

    fx1 = "Wilson-ky-" * string(cont) * ".dat"
    f1 = open(fx1,"w")
    writedlm(f1,[klist re2],"\t")
    close(f1)
end
#------------------------------------------------------------------------------
function main()
    i0 = 1
    for m0 in -1:0.1:1
        Wilsonloop(m0,i0)
        i0 += 1
    end
end
#-----------------------------------------------------------------------------------------------
@time Wilsonloop(1.5,1)
# @time main()
```

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)