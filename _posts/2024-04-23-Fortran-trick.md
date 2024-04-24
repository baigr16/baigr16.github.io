---
title: Fortran踩坑记录(持续更新中)
tags:  Fortran Code 
layout: article
license: true
toc: true
key: a20240423
pageview: true
# cover: /assets/images/Julia/julia-logo.png
header:
  theme: dark
  background: 'linear-gradient(to right, rgb(101, 78, 163), rgb(234, 175, 200))'
article_header:
  type: overlay
  theme: dark
  background_color: false
  background_image: 
    gradient: 'linear-gradient(to right, rgb(101, 78, 163), rgb(234, 175, 200))'
    image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
最近在用`Fortran`写程序的时候，在意想不到的地方踩了坑，这篇Bolg整理一下。
{:.info}
<!--more-->
# 前言
用`Fortran`写程序主打一个方便，只要把公式翻译成代码就行了，但是在码代码的时候还是遇到了一些意想不到的问题，这里整理一下，方便自己查阅闭坑。

- 费米分布函数

在把费米分布转成代码为
```fortran
function fermi(ek)
    ! 费米分布函数
    integer,parameter::dp = kind(1.0)
    real(dp) fermi,ek,kbt
    kbt = 0.01
    fermi = 1.0/(exp(ek/kbt) + 1)  
    ! fermi = ek 
    return
end
```

但这里有$e^x$指数函数，如果$x$足够大的话，如果是在单精度下面(dp = kind(1.0))，此时分布函数就会溢出。不过就算使用双精度(kind(1.0d0))，还是无法完全避免数据溢出这个问题，所以最好的方式就是给定一个阈值，判断在哪些范围内可以使用费米分布函数，其他范围就是0或者1，取决于该量子态是占据还是非占据
```fortran
function fermi(ek)
    !  万万不能直接用费米分布函数，会存在浮点溢出
    integer,parameter::dp = kind(1.0)
    real(dp) fermi,ek,kbt
    kbt = 0.001
    if(ek/kbt>-40 .and. ek/kbt<40) fermi = 1.0/(exp(ek/kbt) + 1) 
    if(ek/kbt <-40) fermi = 1.0
    if(ek/kbt >40) fermi = 0.0
    return
end 
```

- 动态数组赋值

在`Fortran`中经常会用到可变大小的数组，通常会在程序执行过程中来确定数组的大小。但这个时候还是需要对数组进行初始化，但有时候会有错误的初始化，比如

```fortran
real,allocatable::ones(:,:)

allocate(ones(10,10))

do i0 = 1,10
    ones(i0,i0) = 1.0
end do

```
好吧，这是一个错误的示范，因为这里只是对对角元素进行了赋值，但没有对其它位置的元素赋值。所以在后面的程序中如果使用了`ones`这个数组，那么必然就会得到错误的结果，因为其它位置的元素在这里没有赋值，但是`Fortran`对于这种动态类型的数组，不会默认其初始位置都是零，所以那些没有赋值的地方，数据的指向就是奇奇怪怪的地方，计算结果必然是错的。

正确的做法应该是先对确定大小的动态数组的所有元素进行初始化，然后再填入相对应的值
```fortran
real,allocatable::ones(:,:)

allocate(ones(10,10))
ones = 0.0

do i0 = 1,10
    ones(i0,i0) = 1.0
end do

```








# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg){:.border.rounded}{:width="300px" height="300px"}
<div class="card">
  <div class="card__content">
    <div class="card__header">
      <h4>Email</h4>
    </div>
    <p>yxli406@gmail.com</p>
  </div>
</div>