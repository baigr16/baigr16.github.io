---
title: Julia多线程并行加速
tags:  Julia
layout: article
license: true
toc: true
key: a20201223
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
pageview: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
在平时的计算过中,通常会遇到多层的循环嵌套求和,这是一个非常消耗时间的事情,自己现在除了在大型矩阵对角化的时候还是利用Fortran在写程序,平时一些涉及到矢量运算或者动量空间中的一些计算都是在利用Julia做的,那么将这个很耗时的事情变成多线程自然是可以节省很多时间,所以这里就折腾一下如何把一个求和分发到多个线程,然后将所有的结果在放到一起,实现多线程并行,可以节省很多的时间.
{:.info}
<!--more-->
# 简单测试
```julia
using ProgressMeter
using Distributed
@everywhere using SharedArrays, LinearAlgebra
# ----------------------------------------------------------------
@everywhere function fun1(qx::Float64,qy::Float64)::Float64
    kn = 15000
    re::Float64 = 0.0
    for i1 in 0:kn
        for i2 in 0:kn
            kx = i1*pi/kn
            ky = i2*pi/kn
            re = re + kx + ky + qy + qx
        end 
    end
    return re/kn^2
end
#-----------------------------------------
function fun2(kn::Int64)::Float64
    re::Float64 = 0.0
    for i1 in 1:kn
        for i2 in 1:kn
            qx = i1*pi/kn
            qy = i2*pi/kn
            re = re + fun1(qx,qy)
        end
    end
    return re
end
#----------------------------------------
function fun3(kn)
    # 增加8线程，注意别开太多
    addprocs(8 - nprocs());
    println("Running ",nprocs()," processes");
    re = SharedArray(zeros(Float64,1))
    @sync @distributed for i1 in 1:kn
        for i2 in 1:kn
            qx = i1*pi/kn
            qy = i2*pi/kn
            re[1] = re[1] +  fun1(qx,qy)
        end
    end
    return re[1]
end
#-----------------------------------------
function fun4(kn)::SharedArray{Float64,2}
    # 增加8线程，注意别开太多
    addprocs(8 - nprocs());
    println("Running ",nprocs()," processes");
    re = SharedArray(zeros(Float64,kn,kn))
    @sync @distributed for i1 in 1:kn
        for i2 in 1:kn
            qx = i1*pi/kn
            qy = i2*pi/kn
            re[i1,i2] = fun1(qx,qy)
        end
    end
    return re
end
#-------------------------------------------------
function main()
    kn = 15
    @time fun1(0.1,0.1)
    println("Start fun2............")
    @time fun2(kn)
    println("Start fun3............")
    @time fun3(kn)
    println("Start fun4............")
    @time fun4(kn)
end
# ==================================================
main()
```
运行程序的时候需要手动开启多个进程,否则会报错
```shell
julia -p 10 filename.jl
```
如果使用
```shell
julia filename.jl
```
程序则会报错,始终显示函数`fun1`未定义.

# 代码分析
首先将问题描述一下,假设这里要计算一个量,需要进行一些嵌套循环,而这个量在计算的时候又需要其它的参数依赖,这也就是上面的函数`fun1`所描述的问题,他的计算跟变量`qx,qy`是有关系的,而接下来要计算的另外一个量,就是需要对这个`qx,qy`进行撒点,然后再所有的点上对`fun1`进行计算,这也就是下面的三个函数`fun2,fun3,fun4`所描述的计算.

首先来看函数`fun2`
```julia
function fun2(kn::Int64)::Float64
    re::Float64 = 0.0
    for i1 in 1:kn
        for i2 in 1:kn
            qx = i1*pi/kn
            qy = i2*pi/kn
            re = re + fun1(qx,qy)
        end
    end
    return re
end
```
这就是最简单的写法,是个串行执行的程序,只能简单进行循环加和,虽然Julia在这种嵌套循环上已经进行过优化,速度可很客观,但是在现在这个计算情况下他的速度还是很慢的,哪有那么多时间等待程序的执行,可以发现随着`kn`变大,执行时间差不多就是指数上升,所以在这里就需要对这个程序进行优化.

因为我对计算机这方面也算是了解一些并行的概念,首先可以想既然是循环求和,那么可以把这个循环扔到不同的线程上,那么每个线程独立的进行自己的循环,最后将不同线程的结果加和起来就又是想要的结果,所以同时启动多个线程的话,计算速度自然是可以有一个可观的提升,但这里要强调一声,并不是越多的线程就一定能节省更多的时间,首先在并行计算的时候,不同的线程之间也是要进行通讯的,这也是需要时间的,如果你开很多线程进行计算,看起来好像是可以提速了,但是线程之间的通讯又会消耗很多时间,所以可能并不会对速度又多大的提升,所以合理的更具计算量选择要开启的线程数是非常必要的.
{:.warning}

接下来看经过改动的多线程版本`fun3`
```julia
function fun3(kn)::Float64
    # 增加8线程，注意别开太多
    addprocs(8 - nprocs());
    println("Running ",nprocs()," processes");
    re = SharedArray(zeros(Float64,1))
    @sync @distributed for i1 in 1:kn
        for i2 in 1:kn
            qx = i1*pi/kn
            qy = i2*pi/kn
            re[1] = re[1] +  fun1(qx,qy)
        end
    end
    return re[1]
end
```
在这个升级的函数中,最主要的开了多线程来将循环分发到不同的线程上来加速计算.`nprocs()`用来获取当前开启的线程数目,`addprocs()`用来开启一定数目的线程.**这里要说明一下,想要成功的实现多线程并行,一定要安装好程序最开始using的那些库**.这里因为开了多线程的关系,数据的收集方式就要改变一下,所以就有了`re = SharedArray(zeros(Float64,1))`这个语句,因为函数`fun1`的返回值只是一个单纯的数,所以就用这个函数实现一个多线程共享的一维数组,这样可以用这个数组来存储函数`fun1`的结果(我再这里本想找一个共享的变量而不是数组来存储结果,不过我没有找到实现方法,就先用共享数组的方式了).最后利用两个宏`@sync @distributed`来将循环过程分发到不同的线程上,这个操作还是很友好的,不需要再去自己手动设置什么内容.

接下来将`fun3`改动为`fun4`,直接用一个数组来存储所有的结果,因为这个计算结果在平时会用到,其实也就是把`fun3`的共享数组从1维变成不同的大小而已,并没有什么实质性的改动
```julia
function fun4(kn)::SharedArray{Float64,2}
    # 增加8线程，注意别开太多
    addprocs(8 - nprocs());
    println("Running ",nprocs()," processes");
    re = SharedArray(zeros(Float64,kn,kn))
    @sync @distributed for i1 in 1:kn
        for i2 in 1:kn
            qx = i1*pi/kn
            qy = i2*pi/kn
            re[i1,i2] = fun1(qx,qy)
        end
    end
    return re
end
```

在上面并行的版本中,需要对函数`fun1`做一个声明,也就是用宏@everywhere让所有的线程都可以知道这个函数,同样的在程序刚开始,想要使用共享数组的话,同样要对所有的线程声明`@everywhere using SharedArrays, LinearAlgebra`可以让所有的线程都去很好的执行函数并利用共享数组来存储结果.这里还因为知道`fun4`最终结果的具体类型,所以就对函数的最后做了结果类型声明`::SharedArray{Float64,2}`,这样做的好处是可以提升Julia的执行速度.
{.warning}
# 并行结果分析
![png](/assets/images/Julia/julia-mp.png)

上图是执行的结果,从结果中可以看到,对于串行执行的程序(我这里选的参数在我的笔记本上函数`fun1`执行的时间差不多是1s,所以差不多是执行多少次`fun1`那么需要的时间就是多少s)可以发现串行执行的时间确实是最多的,开启了多线程之后,虽然并不是开启几个线程时间就减少几倍,但是时间上差不多是串行的1/6,所以速度上的提升还是很明显的.
