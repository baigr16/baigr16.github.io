---
title: Python简单并行计算(进程并行)
tags: Python Code
layout: article
license: true
toc: true
key: a20200919
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
最近在慢慢复习之前学过的python的内容,毕竟这个语言现成的轮子还是比较多,各种功能实现起来不必自己亲自上手,这里首先学习了一下简单的进程并行计算,但也只是一个很简单的程序.
{:.info}
<!--more-->
# 主程序
```python
# 导入相应的库
import os   
import time
import multiprocessing as pa
from numba import jit # 导入jit用来对循环进行加速
# -------------------------------------------
# 定义一个函数，后面并行运算时需要这个函数
@jit
def sum1(cont):
    c = 0
    for m in range(cont):
        for n in range(cont):
            c = c + m + n
    return c
# --------------------------------------------
# 开启一些进程
for i in range(pa.cpu_count()):  # pa.cpu_count()可以获取当前计算机上最大的进程数
    p = pa.Process(target=sum1,args=(100000*(i + 1),))  # 在这里,循环每进行以一次,就创建一个新的进程
    p.start()

print(pa.cpu_count())
```
> os是一个与系统操作相关的模块,可以对系统的文件以及一系列性质形成读写,包括获取当前的额进程id

> time是一个时间模块,主要用来获取系统时间,在这里则是用来记录程序执行时间

> 从numba中导入jit是用来对python进行加速的,详情可以参考[这里](https://yxli8023.github.io/2020/09/14/Loop-speed.html)


在这个程序中,开启的进程数直接就是计算机上cpu的总数,如果故意将jit加速的装饰去掉,可以直接看到所有正在执行的进程

![png](/assets/images/research/pa1.png)

这里去掉闭包装饰器`@jit`是因为加上它之后,python在计算循环的时候会非常快,从而看不到上面的多进程执行的过程,这里为了展示一下并行的效果,所以注释掉该语句,关于python循环加速可以参考[Julia,Python,Fortran,Mathematica循环计算速度比较](https://yxli8023.github.io/2020/09/14/Loop-speed.html)这篇博客的内容.
{:.success}

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)