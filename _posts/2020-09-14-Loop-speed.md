---
title: Julia,Python,Fortran,Mathematica循环计算速度比较
tags: Julia Code Python Mathematica
layout: article
license: true
toc: true
key: a20200914
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
在平时进行数值计算时，我经常会遇到的两种比较消耗时间计算，其一是对大型矩阵的对角化问题，这是在对哈密顿量厄密矩阵求解本征问题时经常遇到的；另外一个就是循环计算了，这个结构在动量空间计算体系格林函数等其它量时，也会经常遇到。在这里首先对比一下各编程语言对循环的计算速度，然后根据我的经验提供一些加速的方法。
{:.info}
<!--more-->
# 速度对比
在这里只是简单的展示一下计算速度，所以我只是采用了三层循环嵌套，每层的循环数为1000
## python
```python
import time
def f1():
    c = 0
    cont = 1000 #控制循环的次数
    for i in range(cont):
        for j in range(cont):
            for k in range(cont):
                c = c + i + j + k
    return c

t1 = time.time()
print(f1())   # 计算循环求和
t2 = time.time()
print('Timing cost is(s):  ',t2 - t1)
```
> Result is:  1498500000000
> Timing cost is(s):   85.73423385620117

## julia
```julia
function sum1(num::Int128)
    c::Int64 = 0
    for i in 0:num-1
        for j in 0:num-1
            for k in 0:num-1
                c = c + i + j + k
            end
        end
    end
    return c
end
@time sum1(999)
```
> Result is: 1498500000000
> 
> Timing cost is(s): 0.000084 seconds

## Fortran
```fortran
program main
implicit none
integer i,j,k,cont
integer(kind = 8) c
real t1,t2
cont = 999
open(12,file = "time.dat")
call cpu_time(t1)
do i = 0,cont
    do j = 0,cont
        do k = 0,cont
            c = c + i + j + k
        end do
    end do
end do
call cpu_time(t2)
write(12,*)"Result is: ",c
write(12,*)"Timing cost is(s): ",(t2-t1)
close(12)
stop
end program main
```
> Result is:          1498500000000
>
> Timing cost is(s):   0.4992760

## Mathematica
```mathematica
Sum[i + j + k, {i, 0, 999}, {j, 0, 999}, {k, 0, 999}] // AbsoluteTiming
```

![png](/assets/images/research/speed1.png)

**从上面的结果中可以看到，julia在循环计算这方面速度方面的优势还是很明显的。fortran的速度也是很快的，毕竟它本来就是以速度见长的。mathematica的循环计算速度是有些慢的，不过这里强调一下，它是可以进行优化的，而且mathematica的计算优势并不是简单的通过这个循环对比就可以概括的，所以这个对比对mathematica来说并没有太大的意义。python在这里看起来是最慢的，这也符合我对这个语言的认知，它具有它灵活简单的优势，在速度方面就比较慢了，但是后面仍然可以给出python在循环时候的优化方案。**
> julia在语言方面和python具有很大的相似性,而且也可以通过对变量指定特定的类型声明来进行提速,所以我最近几乎都是在学习julia进行编程.但是因为python现在生态比较成熟,可用的轮子比较多,所以如果特别依赖这些轮子的化,可能转julia会有一定的难度.

# 速度优化
## python循环加速
python可以通过编译优化来提高计算速度(我对python的了解并不深入,具体可以参考官网),通过从[numba中导入jit](https://numba.readthedocs.io/en/stable/user/jit.html)来对代码进行优化,官网上说通过这种方式可以实现python的速度和C以及Fortran的速度相比较,下面我就通过这种方式来对python的循环进行速度优化.

这里将循环的次数设置为10000,三层循环之后用来做速度对比在时间上的差异会比较明显
```python
from numba import jit
import time

@jit # 函数闭包
# @jit(nopython=True, parallel=True)
def f1():
    c = 0
    cont = 100000
    for i in range(cont):
        for j in range(cont):
            for k in range(cont):
                c = c + i + j + k
    return c


t1 = time.time()
print(f1())   # 计算循环求和
t2 = time.time()
print('Timing cost is(s):  ',t2 - t1)
```
> 2424547410323587072
> 
> Timing cost is(s):   0.08331894874572754

## julia
```julia
function sum1(num::Int128)
    c::Int64 = 0
    for i in 0:num - 1
        for j in 0:num - 1
            for k in 0:num - 1
                c = c + i + j + k
            end
        end
    end
    return c
end
@time sum1(99999)
```
> 0.794240 seconds
> 2424547410323587072

## fortran
```fortran
program main
implicit none
integer i,j,k,cont
integer(kind = 8) c
real t1,t2
cont = 9999
open(12,file = "time.dat")
call cpu_time(t1)
do i = 0,cont
    do j = 0,cont
        do k = 0,cont
            c = c + i + j + k
        end do
    end do
end do
call cpu_time(t2)
write(12,*)"Result is: ",c
write(12,*)"Timing cost is(s): ",(t2-t1)
close(12)
stop
end program main
```
> Result is:      14998500000000000
>
> Timing cost is(s):    397.8017

**在这里我只简单的对比了julia和python的计算速度,Fortran对循环暂时我还没有找到加速的方法,唯一想到的就是并行,而且fortran随着循环次数的增加,所消耗的时间也是非常长.至于Mathematica来说,想要加速方法还是很多的,比如直接调用一些并行命令,或者利用纯函数进行计算,在这里就不多说了,这里主要是想展示一下过程式编程语言的速度问题,关于mathematica的计算问题,如果有机会我也会单独进行整理.**
# 结语
在这里我只是简单的利用循环,主要对julia和python及fortran的速度进行了对比,其实想想并没有很大的实际价值,不过python的这个速度提升,还是可以参考的.毕竟在实际计算的过程中并不可能单纯的只有简单的循环结构进行数据相加,肯定还是有相应的其它内部操作.我层试过利用Fortran和Julia做与T矩阵杂质散射相关的计算,因为在循环过程中涉及到了矩阵对角化的问题,所以最后在速度上来讲还是Fortran更胜一筹.

最后在这里感谢一下[王泽庆](https://zqw.ink/),这里python的加速建议也是他告诉我的,因为之前我对python的使用率不是很高,所以也并未关注这个加速问题.
