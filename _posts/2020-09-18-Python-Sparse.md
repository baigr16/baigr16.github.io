---
title: Python稀疏矩阵对角化库
tags: Python Code
layout: article
license: true
toc: true
key: a20200918
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
平时在做计算的时候,因为总要对哈密顿量进行对角化,但是这个矩阵一般都是稀疏矩阵,虽然Fortran也有专门的稀疏矩阵对角化的程序,但是使用起来并不是很方便,于是就找到了python专门用来做稀疏矩阵对角化的库,而且利用这个库中的函数,还可以只求特殊要求下的本征值和本征矢量,这也是一个很好的优势,因为通常在做哈密顿量的计算时,也只是关系能量最低的本征态上的问题.
{:.info}
<!--more-->
# 函数介绍
[SciPy](https://docs.scipy.org/doc/scipy/reference/index.html)是一个开源的库,里面包含很多数值计算方面的函数,我也在慢慢学习,这里主要就是介绍一下它对稀疏矩阵对角化的操作.虽然它也是在ARPACK这个Fortran程序之上封装的,但是前面也说过了,如果要使用Fortran还是要花费一定的精力的,毕竟Fortran的接口并不是很好搞定的事情,而python相对来说还是很好用的,所以我也就渐渐转向python的怀抱了,但是不得不说在有些计算方面还是要用Fortran,毕竟计算速度还是杠杠的.

## eigs
```python
scipy.sparse.linalg.eigs(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, OPpart=None)
```

虽然这个函数有许多可选参数,但是提前已经设置好了默认值,所以最简单的使用就是只要提供矩阵$A$就可以得到k个本征值和本征矢量,这里默认k=6.

> A是要求解的矩阵,必须是个方阵,可以是实数型也可以是复数型
> k是一个整数,用来设置你想计算矩阵A的多少个本征值和本征矢量,当然了这个数是不能大于你矩阵的维数的
> sigma可以是实数,也可以是复数,它的设置是为了求解在sigma附近的本征值


> which, [‘LM’ | ‘SM’ | ‘LR’ | ‘SR’ | ‘LI’ | ‘SI’]
> which必须是个字符串,它用来控制你想计算那个位置上的本征值
>‘LM’ : largest magnitude 最大本征值(模)
>‘SM’ : smallest magnitude 最小本征值(模)
>‘LR’ : largest real part 实部最大的本征值
>‘SR’ : smallest real part 实部最小的本征值
>‘LI’ : largest imaginary part 虚部最大的本征值
>‘SI’ : smallest imaginary part 虚部最小的本征值

> maxiter 这个参数是整形的,用来控制计算本征值过程中的最大迭代次数
> tol  浮点型,用来控制计算精度

这个函数可以求解的矩阵是比较一般的,如果所求解的矩阵没有对称性或者其它特殊的性质,通常就使用这个函数来进行对角化.对于没有介绍的参数可以自行参考其网站,因为我再平时的使用中并没有频繁使用,所以就不进行介绍了.
{:.warning}

## eigsh 
```python
scipy.sparse.linalg.eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, mode='normal')
```

**这个函数主要是用来求解实对称矩阵或者是厄密矩阵,我在计算的时候通常是使用这个函数,因为通常我遇到的哈密顿量都是厄密的.**

> sigma 由于这个函数操作的矩阵是有对称性要求的,所以这时候的本征值也有一定的限制,比如厄密矩阵的本征值就一定是实数,所以这时候sigma这个参数只能是个实数.

> which [‘LM’ | ‘SM’ | ‘LA’ | ‘SA’ | ‘BE’] 这个参数设置和eigs是相同的,但不同的是'BE'这个选项使用的条件变了
> If A is a complex hermitian matrix, ‘BE’ is invalid. Which k eigenvectors and eigenvalues to find:
> ‘LM’ : Largest (in magnitude) eigenvalues.
> ‘SM’ : Smallest (in magnitude) eigenvalues.
> ‘LA’ : Largest (algebraic) eigenvalues.
> ‘SA’ : Smallest (algebraic) eigenvalues.
> ‘BE’ : Half (k/2) from each end of the spectrum.

# 程序示例
```python
import numpy as np
from scipy.linalg import eig, eigh
from scipy.sparse.linalg import eigs, eigsh
import time
np.set_printoptions(suppress=True)
# ----------------------------------------------------
cont = 10000
np.random.seed(0)
X = np.random.random((cont,cont)) - 0.5
X = np.dot(X, X.T) #create a symmetric matrix
# -----------------------------------------------------
t1 = time.time()
#evals_all, evecs_all = eigh(X)
#evals_large, evecs_large = eigsh(X, 3, which='LM')
evals_be, evecs_be = eigsh(X, 4, which='BE')
t2 = time.time()
print("Time cost is(s): %s"%(t2 - t1))

```

# 参考
1.[SciPy](https://docs.scipy.org/doc/scipy/reference/index.html)