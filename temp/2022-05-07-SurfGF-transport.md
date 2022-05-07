---
title: 的带方法自能计算
tags: Julia transport 
layout: article
license: true
toc: true
key: a20220505
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
    g2 = a0 * g1 * b0 # 自能
    g3 = im * (g2 - conj(g2)) # 展宽函数
    println(g2)
    return g1,g2,g3
end
#-------------------------------------------------
@time g1,g2,g3 = SurfGF()
```

# 参考
- 1.[输运问题总结与练习(一):二端体系](https://zhuanlan.zhihu.com/p/269595149)