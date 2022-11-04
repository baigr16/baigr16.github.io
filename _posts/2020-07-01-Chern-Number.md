---
title: 两种方法计算Chern Number
tags: Study Julia Python Code Topology
layout: article
license: true
toc: true
pageview: true
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
key: a20200701
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---

计算[Chern数]( https://en.wikipedia.org/wiki/Berry_connection_and_curvature )是最初学习拓扑物理都会遇到的问题，正好在假期空闲的时候自己学习了一下Chern数的数值计算方法，在博客上记录一下希望可以帮助到别人。
{:.info}
<!--more-->
具体的计算方法和细节就不在这里说明了，只要是想学习计算Chern数的肯定了解它在凝聚态物理中的角色，而计算的细节也会在后面的参考文献中给出，只是展示一下结果。
# Julia语言计算Chern number
## Version1
这个方法是直接用定义直接计算的结果，但是可能会遇到波函数规范选择的问题，会导致结果有误，而具体的规范问题，我并不懂，所以一般我会选择第二种方法来计算，也就是后面参考文献中介绍的方法。
```julia
import Pkg
Pkg.add("LinearAlgebra")
Pkg.add("PyPlot")
using LinearAlgebra,PyPlot
function matSet(kx::Float64,ky::Float64)::Matrix{ComplexF64}
    m0::Float64 = -1.0
    t1::Float64 = 1.0
    t2::Float64 = 1.0
    t3::Float64 = 0.5
    # 这里选取的是量子反常Hall效应的模型
    ham = zeros(ComplexF64,2,2)
    ham[1,1] = m0 + 2*t3*sin(kx) + 2*t3*sin(ky) + 2*t2*cos(kx + ky)
    ham[2,2] = -(m0 + 2*t3*sin(kx) + 2*t3*sin(ky) + 2*t2*cos(kx + ky))
    ham[1,2] = 2*t1*cos(kx) - 1im*2*t1*cos(ky)
    ham[2,1] = conj(ham[1,2])  
    return ham
end
#--------------------------------------------------------------------------
function ux(kx::Float64,ky::Float64,ne::Int64)::ComplexF64
    del::Float64 = pi/ne
    #----
    w0 = eigvecs(matSet(kx,ky))[:,1]
    #-----
    wx = eigvecs(matSet(kx + del,ky))[:,1]
    #------
    return  w0'*wx/abs(w0'*wx)
end
#---------------------------------------------------------------------------
function uy(kx::Float64,ky::Float64,ne::Int64)::ComplexF64
    del::Float64 = pi/ne
    #----
    w0 = eigvecs(matSet(kx,ky))[:,1]
    #-----
    wy = eigvecs(matSet(kx,ky + del))[:,1]
    #------
    return  w0'*wy/abs(w0'*wy)
end
#---------------------------------------------------------------------------
function img1(xlist::Array{Float64},ylist::Array{Float64},zlist::Array{ComplexF64})
    zlist = map(imag,zlist)
    #p1 = scatter(xlist,ylist,zlist*20,c=zlist*0.1,edgecolors="b",cmap="Reds")
    p1 = scatter(xlist,ylist,zlist*200,c=zlist,cmap="Reds")
    colorbar(p1)
    xlabel("kx")
    ylabel("ky")
    title("Berry Curvature")
    savefig("Berry Curature.png",bbox_inches="tight",dpi=60)
end
#----------------------------------------------------------------------------
function main(ne::Int64)
    del::Float64 = pi/ne
    kx::Float64 = 0.0
    ky::Float64 = 0.0
    flux::ComplexF64 = 0.0 + 0.0im
    chern_num::ComplexF64 = 0.0 + 0.0im
    kxlist = Float64[]
    kylist = Float64[]
    flist = ComplexF64[]
    for m1 = -ne:ne
        kx = m1*pi/ne
        for m2 = -ne:ne
            append!(kxlist,kx)
            ky = m2*pi/ne
            append!(kylist,ky)
            flux = log((ux(kx,ky,ne)*uy(kx + del,ky,ne))/(ux(kx,ky + del,ne)*uy(kx,ky,ne)))
            append!(flist,flux)
            chern_num = chern_num + flux
        end
    end
    img1(kxlist,kylist,flist)
    return round(real(chern_num/(2.0*pi*1im)))
end
#--------------------------------------------
main(100)
```
## Version2
```julia
using LinearAlgebra,PyPlot
function matSet(kx::Float64,ky::Float64)::Matrix{ComplexF64}
    m0::Float64 = -1.0
    t1::Float64 = 1.0
    t2::Float64 = 1.0
    t3::Float64 = 0.5
    ham = zeros(ComplexF64,2,2)
    ham[1,1] = m0 + 2*t3*sin(kx) + 2*t3*sin(ky) + 2*t2*cos(kx + ky)
    ham[2,2] = -(m0 + 2*t3*sin(kx) + 2*t3*sin(ky) + 2*t2*cos(kx + ky))
    ham[1,2] = 2*t1*cos(kx) - 1im*2*t1*cos(ky)
    ham[2,1] = conj(ham[1,2])  
    return ham
end
#--------------------------------------------------------------------------
function ux(kx::Float64,ky::Float64,ne::Int64)::ComplexF64
    del::Float64 = pi/ne
    #----
    w0 = eigvecs(matSet(kx,ky))[:,1]
    #-----
    wx = eigvecs(matSet(kx + del,ky))[:,1]
    #------
    return  w0'*wx/abs(w0'*wx)
end
#---------------------------------------------------------------------------
function uy(kx::Float64,ky::Float64,ne::Int64)::ComplexF64
    del::Float64 = pi/ne
    #----
    w0 = eigvecs(matSet(kx,ky))[:,1]
    #-----
    wy = eigvecs(matSet(kx,ky + del))[:,1]
    #------
    return  w0'*wy/abs(w0'*wy)
end
#---------------------------------------------------------------------------
function img1(xlist::Array{Float64},ylist::Array{Float64},zlist::Array{ComplexF64})
    zlist = map(imag,zlist)
    #p1 = scatter(xlist,ylist,zlist*20,c=zlist*0.1,edgecolors="b",cmap="Reds")
    p1 = scatter(xlist,ylist,zlist*200,c=zlist,cmap="Reds")
    colorbar(p1)
    xlabel("kx")
    ylabel("ky")
    title("Berry Curvature")
    savefig("Berry Curature.png",bbox_inches="tight",dpi=60)
end
#----------------------------------------------------------------------------
function main(ne::Int64)
    del::Float64 = pi/ne
    kx::Float64 = 0.0
    ky::Float64 = 0.0
    flux::ComplexF64 = 0.0 + 0.0im
    chern_num::ComplexF64 = 0.0 + 0.0im
    kxlist = Float64[]
    kylist = Float64[]
    flist = ComplexF64[]
    for m1 = -ne:ne
        kx = m1*pi/ne
        for m2 = -ne:ne
            append!(kxlist,kx)
            ky = m2*pi/ne
            append!(kylist,ky)
            flux = log((ux(kx,ky,ne)*uy(kx + del,ky,ne))/(ux(kx,ky + del,ne)*uy(kx,ky,ne)))
            append!(flist,flux)
            chern_num = chern_num + flux
        end
    end
    img1(kxlist,kylist,flist)
    return round(real(chern_num/(2.0*pi*1im)))
end
#--------------------------------------------
main(100)
```
# Python Version
在这里要声明，Python版本的是我从[关济寰](http://www.guanjihuan.com/)的网站上复制过来的，必须声明一下版权问题。
```python
import numpy as np
import matplotlib.pyplot as plt
from math import *   # 引入pi, cos等
import time


def hamiltonian(kx, ky):  # 量子反常霍尔QAH模型（该参数对应的陈数为2）
    t1 = 1.0
    t2 = 1.0
    t3 = 0.5
    m = -1.0
    matrix = np.zeros((2, 2))*(1+0j)
    matrix[0, 1] = 2*t1*cos(kx)-1j*2*t1*cos(ky)
    matrix[1, 0] = 2*t1*cos(kx)+1j*2*t1*cos(ky)
    matrix[0, 0] = m+2*t3*sin(kx)+2*t3*sin(ky)+2*t2*cos(kx+ky)
    matrix[1, 1] = -matrix[0,0]
    return matrix


def main():
    start_clock = time.time()
    n = 50  # 积分密度
    delta = 1e-9  # 求导的偏离量
    chern_number = 0  # 陈数初始化
    for kx in np.linspace(-pi, pi, n):
        for ky in np.linspace(-pi, pi, n):
            H = hamiltonian(kx, ky)
            eigenvalue, eigenvector = np.linalg.eig(H)
            vector = eigenvector[:, np.argsort(np.real(eigenvalue))[0]]  # 价带波函数([0]即代表取的是填充的能带波函数)
            # print(np.argsort(np.real(eigenvalue))[0])  # 排序索引（从小到大）
            # print(eigenvalue)  # 排序前的本征值
            # print(np.sort(np.real(eigenvalue)))  # 排序后的本征值（从小到大）
           
            H_delta_kx = hamiltonian(kx+delta, ky) 
            eigenvalue, eigenvector = np.linalg.eig(H_delta_kx)
            vector_delta_kx = eigenvector[:, np.argsort(np.real(eigenvalue))[0]]   # 略偏离kx的波函数

            H_delta_ky = hamiltonian(kx, ky+delta)  
            eigenvalue, eigenvector = np.linalg.eig(H_delta_ky)
            vector_delta_ky = eigenvector[:, np.argsort(np.real(eigenvalue))[0]]  # 略偏离ky的波函数

            H_delta_kx_ky = hamiltonian(kx+delta, ky+delta)  
            eigenvalue, eigenvector = np.linalg.eig(H_delta_kx_ky)
            vector_delta_kx_ky = eigenvector[:, np.argsort(np.real(eigenvalue))[0]]  # 略偏离kx和ky的波函数

            # 价带的波函数的贝里联络(berry connection) # 求导后内积
            A_x = np.dot(vector.transpose().conj(), (vector_delta_kx-vector)/delta)   # 贝里联络Ax（x分量）
            A_y = np.dot(vector.transpose().conj(), (vector_delta_ky-vector)/delta)   # 贝里联络Ay（y分量）
            
            A_x_delta_ky = np.dot(vector_delta_ky.transpose().conj(), (vector_delta_kx_ky-vector_delta_ky)/delta)  # 略偏离ky的贝里联络Ax
            A_y_delta_kx = np.dot(vector_delta_kx.transpose().conj(), (vector_delta_kx_ky-vector_delta_kx)/delta)  # 略偏离kx的贝里联络Ay

            # 贝里曲率(berry curvature)
            F = (A_y_delta_kx-A_y)/delta-(A_x_delta_ky-A_x)/delta

            # 陈数(chern number)
            chern_number = chern_number + F*(2*pi/n)**2
    chern_number = chern_number/(2*pi*1j)
    print('Chern number = ', chern_number)
    end_clock = time.time()
    print('CPU执行时间(min)=', (end_clock-start_clock)/60)

if __name__ == '__main__':
    main()
```
# Julia调用Python计算Chern 数
```julia
#  Package import which required
import Pkg
#Pkg.add("PyPlot")
using LinearAlgebra,PyCall
@pyimport numpy as np
#-------------------------------------------------
function hamiltonian(kx,ky)
    t1 = 1.0
    t2 = 1.0
    t3 = 0.5
    m0 = -1.0
    ham = zeros(ComplexF64,2,2)
    ham[1,1] = m0 + 2*t3*sin(kx) + 2*t3*sin(ky) + 2*t2*cos(kx + ky)
    ham[2,2] = -(m0 + 2*t3*sin(kx) + 2*t3*sin(ky) + 2*t2*cos(kx + ky))
    ham[1,2] = 2*t1*cos(kx) - 1im*2*t1*cos(ky)
    ham[2,1] = 2*t1*cos(kx) + 1im*2*t1*cos(ky)
    return ham
end

function main()
    n = 100  # 积分密度
    delta = 1e-9  # 求导的偏离量
    chern_number = 0  # 陈数初始化
    for kx in np.linspace(-pi, pi, n)
        for ky in np.linspace(-pi, pi, n)
            #-------------------------------------------------
            H = hamiltonian(kx, ky)
            eigenvalue, eigenvector = np.linalg.eig(H)
            vector = eigenvector[:, np.argsort(np.real(eigenvalue))[1]+1]  # 价带波函数([0]即代表取的是填充的能带波函数)
           
            H_delta_kx = hamiltonian(kx + delta, ky) 
            eigenvalue, eigenvector = np.linalg.eig(H_delta_kx)
            vector_delta_kx = eigenvector[:, np.argsort(np.real(eigenvalue))[1]+1]   # 略偏离kx的波函数

            H_delta_ky = hamiltonian(kx, ky + delta)  
            eigenvalue, eigenvector = np.linalg.eig(H_delta_ky)
            vector_delta_ky = eigenvector[:, np.argsort(np.real(eigenvalue))[1]+1]  # 略偏离ky的波函数

            H_delta_kx_ky = hamiltonian(kx + delta, ky + delta)  
            eigenvalue, eigenvector = np.linalg.eig(H_delta_kx_ky)
            vector_delta_kx_ky = eigenvector[:, np.argsort(np.real(eigenvalue))[1]+1]  # 略偏离kx和ky的波函数

            # 价带的波函数的贝里联络(berry connection) # 求导后内积
            A_x = vector'*(vector_delta_kx - vector)/delta   # 贝里联络Ax（x分量）
            A_y = vector'*(vector_delta_ky - vector)/delta   # 贝里联络Ay（y分量）
            
            A_x_delta_ky = vector_delta_ky'*(vector_delta_kx_ky - vector_delta_ky)/delta  # 略偏离ky的贝里联络Ax
            A_y_delta_kx = vector_delta_kx'*(vector_delta_kx_ky - vector_delta_kx)/delta  # 略偏离kx的贝里联络Ay

            # 贝里曲率(berry curvature)
            F = (A_y_delta_kx - A_y)/delta - (A_x_delta_ky - A_x)/delta
            # 陈数(chern number)
            chern_number = chern_number + F*(2*pi/n)^2
        end
    end
    return chern_number = np.round(np.real(chern_number/(2*pi*1im)))
end
@time main()
```

# Julia 并行版
```julia
using DelimitedFiles
using ProgressMeter
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles
# -----------------------------------------------------------------------
@everywhere function matSet(kx::Float64,ky::Float64,m0::Float64)::Matrix{ComplexF64}
    #m0::Float64 = 0.5
    tx::Float64 = 1.0
    ty::Float64 = 1.0
    lamx::Float64 = 1.0
    lamy::Float64 = 1.0
    ham = zeros(ComplexF64,2,2)
    ham[1,1] = (m0 - tx*cos(kx) - ty*cos(ky))
    ham[2,2] = -ham[1,1]
    ham[1,2] = lamx*sin(kx) - im*lamy*sin(ky)
    ham[2,1] = conj(ham[1,2])
    return ham
end
#--------------------------------------------------------------------------
@everywhere function ux(kx::Float64,ky::Float64,ne::Int64,m0::Float64)::ComplexF64
    del::Float64 = pi/ne
    #----
    w0 = eigvecs(matSet(kx,ky,m0))[:,1]
    #-----
    wx = eigvecs(matSet(kx + del,ky,m0))[:,1]
    #------
    return  w0'*wx/abs(w0'*wx)
end
#---------------------------------------------------------------------------
@everywhere function uy(kx::Float64,ky::Float64,ne::Int64,m0::Float64)::ComplexF64
    del::Float64 = pi/ne
    #----
    w0 = eigvecs(matSet(kx,ky,m0))[:,1]
    #-----
    wy = eigvecs(matSet(kx,ky + del,m0))[:,1]
    #------
    return  w0'*wy/abs(w0'*wy)
end
#----------------------------------------------------------------------------
@everywhere function ChernNumber(m0::Float64)
    ne::Int64 = 200
    del::Float64 = pi/ne
    kx::Float64 = 0.0
    ky::Float64 = 0.0
    flux::ComplexF64 = 0.0 + 0.0im
    chern_num::ComplexF64 = 0.0 + 0.0im
    kxlist = Float64[]
    kylist = Float64[]
    flist = ComplexF64[]
    for m1 = 0:2*ne
        kx = m1*pi/ne
        for m2 = 0:2*ne
            append!(kxlist,kx)
            ky = m2*pi/ne
            append!(kylist,ky)
            flux = log((ux(kx,ky,ne,m0)*uy(kx + del,ky,ne,m0))/(ux(kx,ky + del,ne,m0)*uy(kx,ky,ne,m0)))
            append!(flist,flux)
            chern_num = chern_num + flux
        end
    end
    return -round(real(chern_num/(2.0*pi*1im)))
end
#--------------------------------------------
@everywhere function main1()
    ch::Float64 = 0.0
    m0len = 200
    m0list = range(-4,4,length = m0len)
    relist = SharedArray(zeros(Float64,m0len,2))
    @sync @distributed for i0 in 1:m0len
        m0 = m0list[i0]
        ch = ChernNumber(m0)
        relist[i0,1] = m0
        relist[i0,2] = ch
        i0 += 1
    end
    f1 = open("chern-3.dat","w")
    writedlm(f1,relist)
    close(f1)
end
# =================================================
@time main1()

```


# 参考文献
- [Chern Numbers in Discretized Brillouin Zone](https://journals.jps.jp/doi/abs/10.1143/JPSJ.74.1674)
- [Numerical determination of Chern numbers and critical exponents for Anderson localization in tight-binding and related models](https://www.semanticscholar.org/paper/%E2%80%9CNumerical-determination-of-Chern-numbers-and-for-Talkington/77f87a883c13e6047559f9a1132ff7f3dcfab84b)