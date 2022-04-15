---
title: Julia 读取Wannier90hr数据
tags: Julia 
layout: article
license: true
toc: true
key: a20220405
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
最近准备完全使用`Julia`进行科研了，所以就学习一下用该语言进行文件数据的读写如何进行，这里先是一个简单的尝试，就像用`Julia`读取紧束缚的`hr`数据来构建模型，顺便熟悉一下`julia`。
{:.info}
<!--more-->
# 代码展示
```julia
function Wannier90Read(filenam)
    f1 = open(filenam)
    data = readlines(f1)
    Hr = zeros(Float64,length(data[5:end]),7)
    hn = 0
    hopnum = 0
    degen = []
    i0 = 1
    for da in data[2:end] # 直接读取hr的主体数据
        if i0 == 1
            temp = strip(da) # 先去除左右两端空格
            temp = replace(temp,r"\s +"=>",") # 将空格转换为逗号
            temp = split(temp,[',']) # 按照逗号进行分割
            hn = map(x->parse(Float64,x),temp) # 将分割后的数据变化为浮点数
        elseif i0 == 2
            temp = strip(da) # 先去除左右两端空格
            temp = replace(temp,r"\s +"=>",") # 将空格转换为逗号
            temp = split(temp,[',']) # 按照逗号进行分割
            hopnum = map(x->parse(Float64,x),temp) # 将分割后的数据变化为浮点数
        elseif i0==3
            temp = strip(da) # 先去除左右两端空格
            temp = replace(temp,r"\s +"=>",") # 将空格转换为逗号
            temp = split(temp,[',']) # 按照逗号进行分割
            degen = map(x->parse(Float64,x),temp) # 将分割后的数据变化为浮点数
        else
            temp = strip(da) # 先去除左右两端空格
            temp = replace(temp,r"\s +"=>",") # 将空格转换为逗号
            temp = split(temp,[',']) # 按照逗号进行分割
            temp = map(x->parse(Float64,x),temp) # 将分割后的数据变化为浮点数
            Hr[i0 - 3,:] = temp[:]
        end
        i0 += 1
    end
    return hn,hopnum,degen,Hr
end
#-------------------------------------------------------------
function hamset(kx::Float64,ky::Float64,kz::Float64)
    hn,hopnum,degen,Hr = Wannier90Read("wannier90_hr.dat")
    for i0 in 1:hopnum
        H += exp(im*Hr[i0*hn^2,1:3]*[kx,ky,kz])*Hr[i0*hn^2,7:8] # 构建动量空间哈密顿量
    end
end
#-------------------------------------------------------------
a1,a2,a3,a4 = Wannier90Read("wannier90_hr.dat")
hamset(0.1,0.1,0.0)
```