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
    while err > 10^(-9)
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
    return g2
end
#-------------------------------------------------------------------
function conGF()
    # 计算导体区域的格林函数
    H00,H01 = matset()
    H10 = H01'
    hn::Int64 = size(H00)[1] #矩阵维度
    num = hn
    englist = range(-10,10,length = num)
    trans = []
    eye = Diagonal(repeat([1],outer  = hn)) # 构造单位矩阵
    for i0 in 1:num
        eng = englist[i0] + 10^(-6)*1im
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
    println(a1)
end
#--------------------------------------------------------------------
@time main()
```



# 参考
- 1.[输运问题总结与练习(一):二端体系](https://zhuanlan.zhihu.com/p/269595149)