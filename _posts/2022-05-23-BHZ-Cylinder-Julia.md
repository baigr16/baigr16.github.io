---
title: Julia并行计算cylinder边界态
tags: Julia Topology Python
layout: article
license: true
toc: true
key: a20220523
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
之前虽然也写过cylinder计算的代码，但都是利用`Fortran`写的，虽然`Fortran`的计算速度很快，奈何很多简单的操作实现起来实在不太方便，最近干脆全面转`julia`了，虽然速度比不上`Frotran`，但是我可以并行计算呀，`Fortran`的并行没时间，懒的弄了，等有机会再说，这里就用`Julia`并行的计算边界态。
{:.info}
<!--more-->
# 模型
还是用我最熟悉的模型`BHZ+Superconductor`
$$
H(\mathbf{k})=(m_0-t_x\cos k_x-t_y\cos k_y)\sigma_z\tau_+\lambda_x\sin k_x\sigma_xs_z+\lambda_y\sin k_y\sigma_y\tau_z+\Delta(\mathbf{k})s_y\tau_y
$$
具体怎么实现可以查阅我其他的博客，我这里直接就上代码了

# 代码
```julia
@everywhere using SharedArrays, LinearAlgebra,Distributed,DelimitedFiles,Printf
# =================================================
@everywhere function openx(h0::Float64,yn::Int64,ky::Float64)
    hn::Int64 = 8
    # yn::Int64 = 50
    N::Int64 = yn*hn
    m0::Float64 = 1.0
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    txy::Float64 = 2.0
    #-----------------
    dx::Float64 = 0.
    dy::Float64 = -dx
    d0::Float64 = 0.4
    mu::Float64 = 0.0
    dp::Float64 = 0.3
    #h0::Float64 = 0.6 # 层间耦合 
    tp::Float64 = -0. #  inversion breaking
    Ham = zeros(ComplexF64,N,N)
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    g7 = zeros(ComplexF64,hn,hn)
    g1,g2,g3,g4,g5,g6,g7 = gamma()
    for k = 0:yn-1
        if (k == 0)  # Only right block in first line
            for m = 1:hn
                for l = 1:hn
                    Ham[m,l] = (m0-ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l]  - mu*g7[m,l]

                    Ham[m,l + hn] = (-tx*g1[m,l] - im*ax*g2[m,l])/2.0+ dx/2.0*g4[m,l]
                end 
            end 
        elseif ( k==yn-1 ) # Only left block in last line
            for m = 1:hn
                for l = 1:hn
                    Ham[k*hn + m,k*hn + l] = (m0-ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l] - mu*g7[m,l]

                    Ham[k*hn + m,k*hn + l - hn] = -tx*g1[m,l]/2 + im*ax*g2[m,l]/2 + dx/2.0*g4[m,l]
                end
            end
        else
            for m = 1:hn
                for l = 1:hn # k start from 1,matrix block from 2th row
                    Ham[k*hn + m,k*hn + l] = (m0 - ty*cos(ky))*g1[m,l] + ay*sin(ky)*g3[m,l] + (d0 + dy*cos(ky))*g4[m,l] - mu*g7[m,l]

                    Ham[k*hn + m,k*hn + l + hn] = (-tx*g1[m,l] - im*ax*g2[m,l])/2 + dx/2.0*g4[m,l]
                    Ham[k*hn + m,k*hn + l - hn] = -tx*g1[m,l]/2 + im*ax*g2[m,l]/2 + dx/2.0*g4[m,l]
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
    m0::Float64 = 1.0
    tx::Float64 = 2.0
    ty::Float64 = 2.0
    ax::Float64 = 2.0
    ay::Float64 = 2.0
    txy::Float64 = 2.0
    #-----------------
    dx::Float64 = 0.
    dy::Float64 = -dx
    d0::Float64 = 0.4
    dp::Float64 = 0.3
    mu::Float64 = 0.0
#     h0::Float64 = 0.2 # 层间耦合 
    tp::Float64 = -0. #  inversion breaking
    Ham = zeros(ComplexF64,N,N)
    g1 = zeros(ComplexF64,hn,hn)
    g2 = zeros(ComplexF64,hn,hn)
    g3 = zeros(ComplexF64,hn,hn)
    g4 = zeros(ComplexF64,hn,hn)
    g5 = zeros(ComplexF64,hn,hn)
    g6 = zeros(ComplexF64,hn,hn)
    g7 = zeros(ComplexF64,hn,hn)
    g1,g2,g3,g4,g5,g6,g7 = gamma()
    for k = 0:yn-1
        if (k == 0) # Only right block in first line
            for m = 1:hn
                for l = 1:hn
                    Ham[m,l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l] - mu*g7[m,l]

                    Ham[m,l + hn] = (-ty*g1[m,l] - im*ay*g3[m,l])/2 + dy/2.0*g4[m,l]
                end
            end
        elseif ( k==yn-1 ) # Only left block in last line
            for m = 1:hn
                for l = 1:hn
                    Ham[k*hn + m,k*hn + l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l] - mu*g7[m,l]

                    Ham[k*hn + m,k*hn + l - hn] = -ty*g1[m,l]/2 + im*ay*g3[m,l]/2 + dy/2.0*g4[m,l]
                end
            end
        else
            for m = 1:hn
                for l = 1:hn # k start from 1,matrix block from 2th row
                    Ham[k*hn + m,k*hn + l] = (m0-tx*cos(kx))*g1[m,l] + ax*sin(kx)*g2[m,l] + (d0 + dx*cos(kx))*g4[m,l] - mu*g7[m,l]

                    Ham[k*hn + m,k*hn + l + hn] = (-ty*g1[m,l] - im*ay*g3[m,l] )/2 + dy/2.0*g4[m,l]
                    Ham[k*hn + m,k*hn + l - hn] = -ty*g1[m,l]/2 + im*ay*g3[m,l]/2 + dy/2.0*g4[m,l]
                end
            end
        end
    end
    return Ham
end 
#-------------------------------------------------------------------
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
#------------------------------------------------------
@everywhere function cylinder(h0::Float64)
    # h0::Float64 = 0.
    hn::Int64 = 8
    yn::Int64 = 50
    N::Int64 = hn*yn
    ham = zeros(ComplexF64,N,N)
    kn::Int64 = 50
    vals1 = zeros(Float64,2*kn + 1,N)
    vals2 = zeros(Float64,2*kn + 1,N)
    klist = []
    for i1 in -kn:kn
        kx = i1*pi/kn
        append!(klist,kx)
        ham1 = openx(h0,yn,kx)
        ham2 = openy(h0,yn,kx)
        val1 = eigvals(ham1)
        val2 = eigvals(ham2)
        vals1[i1 + kn + 1,:] = map(real,val1[:])
        vals2[i1 + kn + 1,:] = map(real,val2[:])
    end
    fn1 = "ox-" * string(h0) * ".dat"
    fn2 = "oy-" * string(h0) * ".dat"
    f1 = open(fn1,"w")
    f2 = open(fn2,"w")
    klist = (a->(@sprintf "%15.8f" a)).(klist)
    vals1 = (a->(@sprintf "%15.8f" a)).(vals1)
    vals2 = (a->(@sprintf "%15.8f" a)).(vals2)
    writedlm(f1,[klist vals1])
    writedlm(f2,[klist vals2])
    close(f1)
    close(f2)
end
#-------------------------------------------------------
@everywhere function main1()
    @sync @distributed for h0 in -2:0.1:2
        cylinder(h0)
    end
end
#------------------------------------------------------
@time main1()
```

# 绘图
`Julia`画图功能暂时不是很完善，所以就用`Python`来绘图了，下面上绘图代码
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
config = {
"font.size": 30,
"mathtext.fontset":'stix',
"font.serif": ['SimSun'],
}
rcParams.update(config) # Latex 字体设置
#---------------------------------------------------------
def scatterplot1(cont):
    #da1 = "m" + str(cont) + "-pro-ox"  + ".dat"
    #da2 = "m" + str(cont) + "-pro-oy"  + ".dat"
    #da1 = "ox-" + str(cont).rjust(2,'0') + ".dat"
    da1 = "ox-" + str(cont) + ".dat"
    picname = "ox-" + str(cont) + ".png"
    os.chdir(os.getcwd())# 确定用户执行路径
    x0 = []
    y0 = []
    with open(da1) as file:
        da = file.readlines()
        for f1 in da:
            if len(f1) > 3:
                ldos = [float(x) for x in f1.strip().split()]
                x0.append(ldos)
                #y0.append(ldos)
    x0 = np.array(x0)
    plt.figure(figsize=(8,8))
    plt.plot(x0[:,0], x0[:,1:-1], c = 'darkblue', alpha = 0.5)
    plt.plot(x0[:,0], x0[:,int(len(x0[1,:])/2)], c = 'red')
    plt.plot(x0[:,0], x0[:,int(len(x0[1,:])/2) + 1], c = 'red')
    x0min = np.min(x0[:,0])
    x0max = np.max(x0[:,0])
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 30,
             }
    plt.xlim(x0min,x0max)
    plt.ylim(-3,3)
    plt.xlabel(r'$k_y/\pi$',font2)
    plt.ylabel("E",font2)
    #plt.yticks(fontproperties='Times New Roman', size = 15)
    #plt.xticks(fontproperties='Times New Roman', size = 15)
    plt.xticks([-1,0,1],fontproperties='Times New Roman', size = 30)
    plt.yticks([-3,0,3],fontproperties='Times New Roman', size = 30)
    plt.savefig(picname, dpi = 100, bbox_inches = 'tight')
    plt.close()
#---------------------------------------------------------
def scatterplot2(cont):
    #da1 = "m" + str(cont) + "-pro-ox"  + ".dat"
    #da2 = "m" + str(cont) + "-pro-oy"  + ".dat"
    #da1 = "did-oy-" + str(cont).rjust(2,'0') + ".dat"
    da1 = "oy-" + str(cont) + ".dat"
    picname = "oy-" + str(cont) + ".png"
    os.chdir(os.getcwd())# 确定用户执行路径
    x0 = []
    y0 = []
    with open(da1) as file:
        da = file.readlines()
        for f1 in da:
            if len(f1) > 3:
                ldos = [float(x) for x in f1.strip().split()]
                x0.append(ldos)
                #y0.append(ldos)
    x0 = np.array(x0)
    plt.figure(figsize=(8,8))
    plt.plot(x0[:,0], x0[:,1:-1], c = 'darkblue', alpha = 0.5)
    plt.plot(x0[:,0], x0[:,int(len(x0[1,:])/2)], c = 'red')
    plt.plot(x0[:,0], x0[:,int(len(x0[1,:])/2) + 1], c = 'red')
    x0min = np.min(x0[:,0])
    x0max = np.max(x0[:,0])
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 30,
             }
    plt.xlim(x0min,x0max)
    plt.ylim(-3,3)
    plt.xlabel("$k_x/\pi$",font2)
    plt.ylabel("E",font2)
    #plt.yticks(fontproperties='Times New Roman', size = 15)
    #plt.xticks(fontproperties='Times New Roman', size = 15)
    plt.xticks([-1,0,1],fontproperties='Times New Roman', size = 30)
    plt.yticks([-3,0,3],fontproperties='Times New Roman', size = 30)
    plt.savefig(picname, dpi = 100, bbox_inches = 'tight')
    plt.close()
#---------------------------------------------------------
def main():
    for i0 in np.linspace(-2,2,41):
        scatterplot1(format(i0,'.1f'))  
        scatterplot2(format(i0,'.1f')) 
#---------------------------------------------------------
if __name__=="__main__":
    main()
    #scatterplot1(1)
```
这里因为`Julia`在计算过程中数据输出的时候是按照参数的值输出的，所以在绘图的时候需要对脚本做一些小的处理
```python
def main():
    for i0 in np.linspace(-2,2,41):
        scatterplot1(format(i0,'.1f'))  
        scatterplot2(format(i0,'.1f')) 
```
这里将输入的参量进行了格式化，和文件名匹配，再进行绘图。

![png](/assets/images/python/fig5.png)