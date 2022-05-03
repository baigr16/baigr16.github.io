---
title: 格林函数方法学习输运笔记代码整理
tags: transport Julia
layout: article
license: true
toc: true
key: a20220503
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
这里整理一下学习输运时候的代码和笔记。主要的内容都是从[关济寰](https://www.guanjihuan.com/)的网站上学习而来。
{:.info}
<!--more-->
```julia
using LinearAlgebra,DelimitedFiles
#------------------------------------------
function H00(wid::Int64)
    ham = zeros(ComplexF64,wid,wid)
    for i0 in 1:wid - 1 #边界上没有hopping
        ham[i0,i0 + 1] = 1
        ham[i0 + 1,i0] = 1
    end 
    return ham
end
#----------------------------------------
function H01(wid::Int64)
    ham = zeros(ComplexF64,wid,wid)
    for i0 in 1:wid
        ham[i0,i0] = 1
    end 
    return ham
end
#-------------------------------------------
function hamset(wid::Int64,len::Int64) # 定义实空间完整的哈密顿量
    ham = zeros(ComplexF64,wid * len,wid * len)
    for ix in 1:len
        for iy in 1:wid - 1
            ham[(ix - 1) * wid + iy, (ix - 1) * wid + iy + 1] = 1 # 沿着len方向上的hopping
            ham[(ix - 1) * wid + iy + 1, (ix - 1) * wid + iy] = 1 # 共轭操作
        end
    end
    for ix in 1:len - 1
        for iy in 1:wid
            ham[(ix - 1) * wid + iy, ix * wid + iy] = 1 # 沿着wid方向上的hopping,这里在边界上自然是不存在hopping的
            ham[ix * wid + iy, (ix - 1) * wid + iy] = 1 # 共轭操作
        end
    end
    return ham
end 
#---------------------------------------------------
function diagway(fermi::ComplexF64,ham) 
# 通过直接求逆的方式来得到哈密顿量对应的个人林函数
hn = size(ham)[1] # 获取矩阵的维度,因为这里是方阵所以就取了第一个量
eye = zeros(ComplexF64,hn,hn) # 构建单位矩阵
for i0 in 1:hn
    eye[i0,i0] = 1
end
green = inv(fermi * eye - ham)
return green
end
#---------------------------------------------------
function Dyson(fermi::ComplexF64,h00,h01,len1)
    hn = size(h00)[1]
    eye = zeros(Float64,hn,hn) # 构建单位矩阵
    gn = zeros(Float64,hn,hn) 
    gn0 = zeros(Float64,hn,hn) 
    for i0 in 1:hn
        eye[i0,i0] = 1
    end
    for i0 in 1:len1 # 开始沿着len1的方向利用Dyson方程求解格林函数
        if i0 == 1 # 最左边，这个位置我们通常考虑加入电极，这个时候就可以解的有电极存在时的格林函数
            gn = inv(fermi * eye - h00) # 这里先没有考虑电极，所以就没有减去左端电极的自能
            gn0 = gn
        elseif i0 != len1 # 此时所以导体内部
            gn = inv(fermi * eye  - h00 - transpose(h01) * gn * h01)
            gn0 = gn0 * h01 * gn
        else  # 这是导体右边界位置
            gn = inv(fermi * eye - h00 - transpose(h01) * gn * h01)
            gn0 = gn0 * h01 * gn
        end
    end 
    return gn0
end
#------------------------------------------------------------
function main()
    wid = 4
    len1 = 200
    h00 = H00(wid)
    h01 = H01(wid)
    ham = hamset(wid,len1)
    fermi = 0.1 + 0.01 * 1im # 计算格林函数加上小的虚部

    @time g1 = diagway(fermi,ham)
    @time g2 = Dyson(fermi,h00,h01,len1)
    println(g1[1:wid,(len1 -1) * wid + 1:(len1 - 1) * wid + wid])
    println(g2)
end
#----------------------------------------------------
@time main()
```

# 电导计算
```julia
using LinearAlgebra,DelimitedFiles
#------------------------------------------
function H00(wid::Int64)
    ham = zeros(ComplexF64,wid,wid)
    for i0 in 1:wid - 1 #边界上没有hopping
        ham[i0,i0 + 1] = 1
        ham[i0 + 1,i0] = 1
    end 
    return ham
end
#----------------------------------------
function H01(wid::Int64)
    ham = zeros(ComplexF64,wid,wid)
    for i0 in 1:wid
        ham[i0,i0] = 1
    end 
    return ham
end
#-----------------------------------------
function transmatrix(fermi,h00,h01,dim) # 这里dim表示子块矩阵h00的维度
    eye = zeros(Float64,dim,dim)
    trans = zeros(ComplexF64,2 * dim, 2 * dim)
    for i0 in 1:dim
        eye[i0,i0] = 1
    end
    trans[1:dim,1:dim] = (inv(h01) * (fermi * eye - h00))[:,:]
    trans[1:dim,dim + 1:2 * dim] = (-inv(h01) * conj(transpose(h01)))[:,:]
    trans[dim + 1:2 * dim,1:dim] = eye[:,:]
    trans[dim + 1:2 * dim,dim + 1:2 * dim] = zeros(Float64,dim,dim)[:,:]
    return trans    
end
#-------------------------------------------------------
function gflead(fermi,h00,h01,dim)
    eye = zeros(Float64,dim,dim)
    for i0 in 1:dim
        eye[i0,i0] = 1
    end
    trans = transmatrix(fermi,h00,h01,dim)
    val,vec = eig(trans)
    temp = zeros(ComplexF64,2 * dim, 2 * dim)
    for i0 in 1:dim
        temp[:,i0] = vec[:,i0]
        i0 += 1
    end
    s1 = temp[dim:2 * dim, 1:dim]
    s2 = temp[1:dim, 1:dim]
    s3 = temp[dim:2 * dim, dim:2 * dim]
    s4 = temp[1:dim, dim:2 * dim]
    rlsurf = inv(fermi * eye - h00 - h01 * s2 * inv(s1))
    llsurf = inv(fermi * eye - h00 - transpose(conj(h01)) * s3 * inv(s4))
    return rlsurf,llsurf
end
#------------------------------------------------------
function selfE(fermi,h00,h01,dim)
    rlsurf,llsurf = gflead(fermi,h00,h01,dim)
    rlse = h01 * rlsurf * transpose(conj(h01))
    llse = transpose(conj(h01)) * llsurf * h01
    return rlse,llse
end
#-------------------------------------------------------
function conductance(fermi,h00,h01,nx = 300) # 计算电导
    dim = size(h00)[1]
    eye = zeros(Float64,dim,dim)
    for i0 in 1:dim
        eye[i0,i0] = 1
    end
    rlse,llse = selfE(fermi,h00,h01,dim)
    for ix in 1:nx
        if ix == 1
            gn = inv(fermi * eye)
            gn0 = gn
        elseif ix != nx
            gn = inv(fermi * eye - h00 - transpose(conj(h01)) * gn * h01)
            gn0 = gn0 * gn
        else
            gn = inv(fermi * eye - h00 - transpose(conj(h01)) * gn * h01 - rlse)
            gn0 = gn0 * h01 * gnn
        end
    end
    gamr = (rlse - transpose(conj(rlse))) * 1im
    gaml = (llse - transpose(conj(llse))) * 1im 
    trans = tr(gaml * gn0 * gamr * transpose(conj(gn0)))
    return trans
end
#---------------------------------------------------
function plotcon(fermiarray,h00,h01)
    dim = size(fermiarray)[1]
    cond = zeros(Float64,dim)
    i0 = 1
    for fermi in fermiarray
        cond0 = real(conductance(fermi + 0.001im,h00,h01))
        cond[i0] = cond0
        i0 += 1
    end
    plot(fermiarray,cond)
end
#-------------------------------------------------
function main()
    wid = 5
    h00 = H00(wid)
    h01 = H01(wid)
    fermiarray = range(-4,4,length = 100)
    # plotcon(fermiarray,h00,h01)
    transmatrix(0.1,h00,h01,10)
end
#-----------------------------------------
@time main()
```
**计算代码有误**

