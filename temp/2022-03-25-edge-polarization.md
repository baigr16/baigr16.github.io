---
title: 边界极化计算
tags: Topology Julia 
layout: article
license: true
toc: true
key: a20220325
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

{:.info}
<!--more-->
# 理论解释

# 代码
## 串行
```julia
# using LinearAlgebra,PyPlot,DelimitedFiles
# s_x\sigma_0\tau_z
using LinearAlgebra,DelimitedFiles
# ==================================================
function Wannier(h0::Float64,yn::Int64,dir::Int64)
    hn::Int64 = 8
#     yn::Int64 = 50
    N::Int64 = yn*hn
    nx::Int64 = 200   # 在离散方向上的撒点数目
    Kx = range(0,2,length = nx)
    ham = zeros(ComplexF64,N,N)
    Noccu::Int64 = Int(N/2) # 占据态数目为哈密顿量矩阵维度一半
    Wave = zeros(ComplexF64,N,N,nx) # 存储每个离散点上的本征矢量（波函数）
    Wan = zeros(ComplexF64,Noccu,Noccu)# 存储Wannier占据态波函数，占据态数目是哈密顿量矩阵的一半
    ang = zeros(Float64,Noccu) # 对每个离散的点，都会对应一组占据态的Wannier本征值  
    F = zeros(ComplexF64,Noccu,Noccu)
    for ix in 1:nx    # 对离散的ki进行Wilson loop构建
        kx = Kx[ix]*pi
        if dir == 0
            ham = openx(h0,yn,kx)
        else
            ham = openy(h0,yn,kx)
        end
        val,vec = eigen(ham)
        Wave[:,:,ix] = vec[:,:]  # 首先存储所有点上的本征矢量
    end
    Wave[:,:,nx] = Wave[:,:,1]  # 首尾相连

    for i1 in 1:Noccu
        F[i1,i1] = 1
    end
    
    for i1 in 1:nx-1
        for i2 in 1:Noccu
            for i3 in 1:Noccu
                Wan[i2,i3] = Wave[:,i2,i1]'*Wave[:,i3,i1 + 1]   # 计算波函数交叠
            end
        end
        F = F*Wan
    end
    F = -im*log(F)
    return F # 得到Wannier哈密顿量
end
# =================================================
function openx(h0::Float64,yn::Int64,ky::Float64)
    hn::Int64 = 8
    # yn::Int64 = 50
    N::Int64 = yn*hn
    m0::Float64 = 1.
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    #-----------------
    dx::Float64 = 0.
    dy::Float64 = -dx
    d0::Float64 = 0.4
    #h0::Float64 = 0.6 # 层间耦合 
    tp::Float64 = -0. #  inversion breaking
    Ham = zeros(ComplexF64,N,N)
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    g1,g2,g3,g4,g5,g6 = Pauli()
    for k = 0:yn-1
        if (k == 0)  # Only right block in first line
            for m = 1:hn
                for l = 1:hn
                    Ham[m,l] = (m0-ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l] + h0*g5[m,l] + tp*g6[m,l]

                    Ham[m,l + hn] = (-tx*g1[m,l] - im*ax*g2[m,l])/2.0+ dx/2.0*g4[m,l]
                end 
            end 
        elseif ( k==yn-1 ) # Only left block in last line
            for m = 1:hn
                for l = 1:hn
                    Ham[k*hn + m,k*hn + l] = (m0-ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l - hn] = -tx*g1[m,l]/2 + im*ax*g2[m,l]/2 + dx/2.0*g4[m,l]
                end
            end
        else
            for m = 1:hn
                for l = 1:hn # k start from 1,matrix block from 2th row
                    Ham[k*hn + m,k*hn + l] = (m0 - ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l + hn] = (-tx*g1[m,l] - im*ax*g2[m,l])/2 + dx/2.0*g4[m,l]
                    Ham[k*hn + m,k*hn + l - hn] = -tx*g1[m,l]/2 + im*ax*g2[m,l]/2 + dx/2.0*g4[m,l]
                end
            end
        end
    end
    return Ham
    end 
# ==========================================================
function openy(h0::Float64,yn::Int64,kx::Float64)
    hn::Int64 = 8
    # yn::Int64 = 50
    N::Int64 = yn*hn
    m0::Float64 = 1.
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    #-----------------
    dx::Float64 = 0.
    dy::Float64 = -dx
    d0::Float64 = 0.4
#     h0::Float64 = 0.2 # 层间耦合 
    tp::Float64 = -0. #  inversion breaking
    Ham = zeros(ComplexF64,N,N)
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    g1,g2,g3,g4,g5,g6 = Pauli()
    for k = 0:yn-1
        if (k == 0) # Only right block in first line
            for m = 1:hn
                for l = 1:hn
                    Ham[m,l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[m,l + hn] = (-ty*g1[m,l] - im*ay*g3[m,l])/2 + dy/2.0*g4[m,l]
                end
            end
        elseif ( k==yn-1 ) # Only left block in last line
            for m = 1:hn
                for l = 1:hn
                    Ham[k*hn + m,k*hn + l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l - hn] = -ty*g1[m,l]/2 + im*ay*g3[m,l]/2 + dy/2.0*g4[m,l]
                end
            end
        else
            for m = 1:hn
                for l = 1:hn # k start from 1,matrix block from 2th row
                    Ham[k*hn + m,k*hn + l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l + hn] = (-ty*g1[m,l] - im*ay*g3[m,l] )/2 + dy/2.0*g4[m,l]
                    Ham[k*hn + m,k*hn + l - hn] = -ty*g1[m,l]/2 + im*ay*g3[m,l]/2 + dy/2.0*g4[m,l]
                end
            end
        end
    end
    return Ham
    end 
# ===================================================
function Pauli()
    hn::Int64 = 8
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    #---------Kinetic energy
    g1[1,1] = 1
    g1[2,2] = -1
    g1[3,3] = 1
    g1[4,4] = -1
    g1[5,5] = -1
    g1[6,6] = 1
    g1[7,7] = -1
    g1[8,8] = 1
    #----------SOC-x
    g2[1,2] = 1
    g2[2,1] = 1
    g2[3,4] = -1
    g2[4,3] = -1
    g2[5,6] = 1
    g2[6,5] = 1
    g2[7,8] = -1
    g2[8,7] = -1
    #---------SOC-y
    g3[1,2] = -im
    g3[2,1] = im
    g3[3,4] = -im
    g3[4,3] = im
    g3[5,6] = im 
    g3[6,5] = -im
    g3[7,8] = im
    g3[8,7] = -im
    #-------------------SC pairing
    g4[1,7] = -1
    g4[2,8] = -1
    g4[3,5] = 1
    g4[4,6] = 1
    g4[7,1] = -1
    g4[8,2] = -1
    g4[5,3] = 1
    g4[6,4] = 1
    #-------------------- layer couple
    g5[1,3] = 1
    g5[2,4] = 1
    g5[3,1] = 1
    g5[4,2] = 1
    g5[5,7] = -1
    g5[6,8] = -1
    g5[7,5] = -1
    g5[8,6] = -1
    #--------------------
    g6[1,2] = 1
    g6[2,1] = 1
    g6[3,4] = 1
    g6[4,3] = 1
    g6[5,6] = -1
    g6[6,5] = -1
    g6[7,8] = -1
    g6[8,7] = -1
    return g1,g2,g3,g4,g5,g6
    end
# ==================================================
function Wannier_spectrum()
# 开边界方向Wannier谱的计算
    h0n::Int64 = 10
    yn::Int64 = 10 # 开边界方向原胞数
    hn::Int64 = 8
    N::Int64 = Int(yn*hn/2)  # 占据态数目
    vals = zeros(Float64,h0n + 1,N)
    hvals = []
    f2 = open("process.dat","w")
    for i1 = 0:h0n
        writedlm(f2,i1/h0n)
        vals[i1 + 1,:] = Wilson(i1/h0n,yn)
        append!(hvals,i1/h0n)
    end
    f1 = open("test2.dat","w")
    writedlm(f1,[hvals vals])
    close(f1)
    close(f2)
    return hvals,vals
end
# ==================================================
function rho(yn::Int64,h0::Float64,dir::Int64)
    # 计算所有格点上的极化分布
    # yn 开边界方向原胞数. h0 层间耦合强度
    hn::Int64 = 8
    N::Int64 = yn*hn
    Nocc::Int64 = Int(N/2)
    kn::Int64 = 50
    re1 = zeros(Float64,yn)
    ham_wan = Wannier(h0,yn,dir) # 得到Wannier哈密顿量
    val2,vec2 = eigen(ham_wan)
    val2 = map(real,val2)/(2*pi)
    for ik = -kn:kn
        kx = ik*pi/kn
        if dir == 0
            println("openx")
            ham = openx(h0,yn,kx)
        else
            println("openy")
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
# =================================================
function edge_pro()
    # 边界极化
    yn::Int64 = 40     # 开边界方向原胞数
    h0n::Int64 = 50
    rho_all = zeros(Float64,h0n,yn)
    h0list = []
    re1::Float64 = 0.0
    for ih = 1:h0n
        rho_all[ih,:] = rho(yn,ih/h0n,0) # 计算出每个格点上的极化强度
        append!(h0list,ih/h0n)
    end

    f1 = open("m1-pro-all-ox.dat","w")
    writedlm(f1,[h0list rho_all])
    close(f1)

    f2 = open("m1-pro-ox.dat","w")
    for ih = 1:h0n
        re1 = 0
        for iy = 1:Int(yn/2)
            re1 += rho_all[ih,iy]  # 计算边界上的极化强度大小,对一半边界格点求和
        end
        re1 = map(x->Int(round(x*2))/2,re1)
        writedlm(f2,[h0list[ih] re1])
    end
    close(f2)
    #--------------------------------------------------
    h0list = []
    for ih = 1:h0n
        rho_all[ih,:] = rho(yn,ih/h0n,2) # 计算出每个格点上的极化强度
        append!(h0list,ih/h0n)
    end

    f1 = open("m1-pro-all-oy.dat","w")
    writedlm(f1,[h0list rho_all])
    close(f1)

    f2 = open("m1-pro-oy.dat","w")
    for ih = 1:h0n
        re1 = 0
        for iy = 1:Int(yn/2)
            re1 += rho_all[ih,iy]  # 计算边界上的极化强度大小,对一半边界格点求和
        end
        re1 = map(x->Int(round(x*2))/2,re1)
        writedlm(f2,[h0list[ih] re1])
    end
    close(f2)
end
# ====================================================
function main1()
    # 计算不同格点位置处的rho
    yn::Int64 = 40  # 开边界方向原胞数
    h0::Float64 = 0.6 # 层间耦合大小
    re1::Float64 = 0
    iylist = []
    relist = []
    for iy = 1:yn
        re1 = rho(yn,h0,iy)
        append!(iylist,iy)
        append!(relist,re1)
    end
    f1 = open("edge_rho.dat","w")
    writedlm(f1,[iylist relist/max(relist)]) # 对结果进行缩放
    close(f1)
end
# ====================================================
@timed edge_pro()
```
## 并行
```julia
# using LinearAlgebra,PyPlot,DelimitedFiles
# s_x\sigma_0\tau_z
using DelimitedFiles
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles
# ==================================================
@everywhere function Wannier(h0::Float64,yn::Int64,dir::Int64)
    hn::Int64 = 8
#     yn::Int64 = 50
    N::Int64 = yn*hn
    nx::Int64 = 200   # 在离散方向上的撒点数目
    Kx = range(0,2,length = nx)
    ham = zeros(ComplexF64,N,N)
    Noccu::Int64 = Int(N/2) # 占据态数目为哈密顿量矩阵维度一半
    Wave = zeros(ComplexF64,N,N,nx) # 存储每个离散点上的本征矢量（波函数）
    Wan = zeros(ComplexF64,Noccu,Noccu)# 存储Wannier占据态波函数，占据态数目是哈密顿量矩阵的一半
    ang = zeros(Float64,Noccu) # 对每个离散的点，都会对应一组占据态的Wannier本征值  
    F = zeros(ComplexF64,Noccu,Noccu)
    for ix in 1:nx    # 对离散的ki进行Wilson loop构建
        kx = Kx[ix]*pi
        if dir == 0
            ham = openx(h0,yn,kx)
        else
            ham = openy(h0,yn,kx)
        end
        val,vec = eigen(ham)
        Wave[:,:,ix] = vec[:,:]  # 首先存储所有点上的本征矢量
    end
    Wave[:,:,nx] = Wave[:,:,1]  # 首尾相连

    for i1 in 1:Noccu
        F[i1,i1] = 1
    end
    
    for i1 in 1:nx-1
        for i2 in 1:Noccu
            for i3 in 1:Noccu
                Wan[i2,i3] = Wave[:,i2,i1]'*Wave[:,i3,i1 + 1]   # 计算波函数交叠
            end
        end
        F = F*Wan
    end
    F = -im*log(F)
    return F # 得到Wannier哈密顿量
end
# =================================================
@everywhere function openx(h0::Float64,yn::Int64,ky::Float64)
    hn::Int64 = 8
    # yn::Int64 = 50
    N::Int64 = yn*hn
    m0::Float64 = -1.
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    h4::Float64 = 2.0
    #-----------------
    dx::Float64 = 0.
    dy::Float64 = -dx
    d0::Float64 = 0.4
    #h0::Float64 = 0.6 # 层间耦合 
    tp::Float64 = -0. #  inversion breaking
    Ham = zeros(ComplexF64,N,N)
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    g1,g2,g3,g4,g5,g6 = Pauli()
    for k = 0:yn-1
        if (k == 0)  # Only right block in first line
            for m = 1:hn
                for l = 1:hn
                    Ham[m,l] = (m0-ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l] + h0*g5[m,l] + tp*g6[m,l]

                    Ham[m,l + hn] = (-tx*g1[m,l] - im*ax*g2[m,l])/2.0+ dx/2.0*g4[m,l] - h4*cos(ky)*g1[m,l]
                end 
            end 
        elseif ( k==yn-1 ) # Only left block in last line
            for m = 1:hn
                for l = 1:hn
                    Ham[k*hn + m,k*hn + l] = (m0-ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l - hn] = -tx*g1[m,l]/2 + im*ax*g2[m,l]/2 + dx/2.0*g4[m,l] - h4*cos(ky)*g1[m,l]
                end
            end
        else
            for m = 1:hn
                for l = 1:hn # k start from 1,matrix block from 2th row
                    Ham[k*hn + m,k*hn + l] = (m0 - ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l + hn] = (-tx*g1[m,l] - im*ax*g2[m,l])/2 + dx/2.0*g4[m,l] - h4*cos(ky)*g1[m,l]
                    Ham[k*hn + m,k*hn + l - hn] = -tx*g1[m,l]/2 + im*ax*g2[m,l]/2 + dx/2.0*g4[m,l] - h4*cos(ky)*g1[m,l]
                end
            end
        end
    end
    return Ham
    end 
# ==========================================================
@everywhere function openy(h0::Float64,yn::Int64,kx::Float64)
    hn::Int64 = 8
    # yn::Int64 = 50
    N::Int64 = yn*hn
    m0::Float64 = -1.
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    h4::Float64 = 2.0
    #-----------------
    dx::Float64 = 0.
    dy::Float64 = -dx
    d0::Float64 = 0.4
#     h0::Float64 = 0.2 # 层间耦合 
    tp::Float64 = -0. #  inversion breaking
    Ham = zeros(ComplexF64,N,N)
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    g1,g2,g3,g4,g5,g6 = Pauli()
    for k = 0:yn-1
        if (k == 0) # Only right block in first line
            for m = 1:hn
                for l = 1:hn
                    Ham[m,l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[m,l + hn] = (-ty*g1[m,l] - im*ay*g3[m,l])/2 + dy/2.0*g4[m,l] - h4*cos(kx)*g1[m,l]
                end
            end
        elseif ( k==yn-1 ) # Only left block in last line
            for m = 1:hn
                for l = 1:hn
                    Ham[k*hn + m,k*hn + l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l - hn] = -ty*g1[m,l]/2 + im*ay*g3[m,l]/2 + dy/2.0*g4[m,l] - h4*cos(kx)*g1[m,l]
                end
            end
        else
            for m = 1:hn
                for l = 1:hn # k start from 1,matrix block from 2th row
                    Ham[k*hn + m,k*hn + l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l]+ h0*g5[m,l] + tp*g6[m,l]

                    Ham[k*hn + m,k*hn + l + hn] = (-ty*g1[m,l] - im*ay*g3[m,l] )/2 + dy/2.0*g4[m,l] - h4*cos(kx)*g1[m,l]
                    Ham[k*hn + m,k*hn + l - hn] = -ty*g1[m,l]/2 + im*ay*g3[m,l]/2 + dy/2.0*g4[m,l] - h4*cos(kx)*g1[m,l]
                end
            end
        end
    end
    return Ham
    end 
# ===================================================
@everywhere function Pauli()
    hn::Int64 = 8
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    #---------Kinetic energy
    g1[1,1] = 1
    g1[2,2] = -1
    g1[3,3] = 1
    g1[4,4] = -1
    g1[5,5] = -1
    g1[6,6] = 1
    g1[7,7] = -1
    g1[8,8] = 1
    #----------SOC-x
    g2[1,2] = 1
    g2[2,1] = 1
    g2[3,4] = -1
    g2[4,3] = -1
    g2[5,6] = 1
    g2[6,5] = 1
    g2[7,8] = -1
    g2[8,7] = -1
    #---------SOC-y
    g3[1,2] = -im
    g3[2,1] = im
    g3[3,4] = -im
    g3[4,3] = im
    g3[5,6] = im 
    g3[6,5] = -im
    g3[7,8] = im
    g3[8,7] = -im
    #-------------------SC pairing
    g4[1,7] = -1
    g4[2,8] = -1
    g4[3,5] = 1
    g4[4,6] = 1
    g4[7,1] = -1
    g4[8,2] = -1
    g4[5,3] = 1
    g4[6,4] = 1
    #-------------------- layer couple
    g5[1,3] = 1
    g5[2,4] = 1
    g5[3,1] = 1
    g5[4,2] = 1
    g5[5,7] = -1
    g5[6,8] = -1
    g5[7,5] = -1
    g5[8,6] = -1
    #--------------------
    g6[1,2] = 1
    g6[2,1] = 1
    g6[3,4] = 1
    g6[4,3] = 1
    g6[5,6] = -1
    g6[6,5] = -1
    g6[7,8] = -1
    g6[8,7] = -1
    return g1,g2,g3,g4,g5,g6
    end
# ==================================================
@everywhere function Wannier_spectrum()
# 开边界方向Wannier谱的计算
    h0n::Int64 = 10
    yn::Int64 = 10 # 开边界方向原胞数
    hn::Int64 = 8
    N::Int64 = Int(yn*hn/2)  # 占据态数目
    vals = zeros(Float64,h0n + 1,N)
    hvals = []
    f2 = open("process.dat","w")
    for i1 = 0:h0n
        writedlm(f2,i1/h0n)
        vals[i1 + 1,:] = Wilson(i1/h0n,yn)
        append!(hvals,i1/h0n)
    end
    f1 = open("test2.dat","w")
    writedlm(f1,[hvals vals])
    close(f1)
    close(f2)
    return hvals,vals
end
# ==================================================
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
# =================================================
@everywhere function edge_pro()
    # 边界极化
    yn::Int64 = 40     # 开边界方向原胞数
    h0n::Int64 = 50
    # rho_all = zeros(Float64,h0n,yn)
    rho_all = SharedArray(zeros(Float64,h0n,yn))
    h0list = SharedArray(zeros(Float64,h0n,1))
    re1::Float64 = 0.0
    @sync @distributed for ih = 1:h0n
        rho_all[ih,:] = rho(yn,ih/h0n,0) # 计算出每个格点上的极化强度
        h0list[ih,1] = ih/h0n
        #append!(h0list,ih/h0n)
    end

    f1 = open("m1-pro-all-ox.dat","w")
    writedlm(f1,[h0list rho_all])
    close(f1)

    f2 = open("m1-pro-ox.dat","w")
    for ih = 1:h0n
        re1 = 0
        for iy = 1:Int(yn/2)
            re1 += rho_all[ih,iy]  # 计算边界上的极化强度大小,对一半边界格点求和
        end
        re1 = map(x->Int(round(x*2))/2,re1)
        writedlm(f2,[h0list[ih] re1])
    end
    close(f2)
    #--------------------------------------------------
    h0list = SharedArray(zeros(Float64,h0n,1))
    @sync @distributed for ih = 1:h0n
        rho_all[ih,:] = rho(yn,ih/h0n,2) # 计算出每个格点上的极化强度
        h0list[ih,1] = ih/h0n
        #append!(h0list,ih/h0n)
    end

    f1 = open("m1-pro-all-oy.dat","w")
    writedlm(f1,[h0list rho_all])
    close(f1)

    f2 = open("m1-pro-oy.dat","w")
    for ih = 1:h0n
        re1 = 0
        for iy = 1:Int(yn/2)
            re1 += rho_all[ih,iy]  # 计算边界上的极化强度大小,对一半边界格点求和
        end
        re1 = map(x->Int(round(x*2))/2,re1)
        writedlm(f2,[h0list[ih] re1])
    end
    close(f2)
end
# ====================================================
@everywhere function main1()
    # 计算不同格点位置处的rho
    yn::Int64 = 40  # 开边界方向原胞数
    h0::Float64 = 0.6 # 层间耦合大小
    re1::Float64 = 0
    iylist = []
    relist = []
    for iy = 1:yn
        re1 = rho(yn,h0,iy)
        append!(iylist,iy)
        append!(relist,re1)
    end
    f1 = open("edge_rho.dat","w")
    writedlm(f1,[iylist relist/max(relist)]) # 对结果进行缩放
    close(f1)
end
# ====================================================
@time edge_pro()
```
# 参考


