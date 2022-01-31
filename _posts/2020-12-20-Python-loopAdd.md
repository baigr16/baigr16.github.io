---
title: Python 循环加速实例
tags: Python 
layout: article
license: true
toc: true
key: a20201220
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
之前在[Julia,Python,Fortran,Mathematica循环计算速度比较](https://yxli8023.github.io/2020/09/14/Loop-speed.html)博客中，我简单的对集中编程语言的循环进行了对比，虽然没什么太大的使用价值，不过对我对自己写代码时候的速度考虑开了一个不错的头，所以这里就把自己利用循环加速改写了一个实例展示出来，看看到底效率如何。
{:.info}
<!--more-->
# 循环加速并行矩阵运算
在前面的博客中，我只是对多重循环中进行了简单的求和，但是在实际的运算中可定不会那么的简单，所以我首先想测试的是如果循环里面是矩阵的运算，速度到底会不会提升的很好。
```python
from numba import jit
import time
import numpy as np

@jit # 函数闭包
#@jit(nopython=True, parallel=True)
def f1():
    c = 0
    cont = 100000
    for i in range(cont):
        for j in range(cont):
            for k in range(cont):
                c = c + i + j + k
    return c

# @jit
#@jit(nopython=True, parallel=True)
@jit(parallel=True)
def f2():
    ndim = 100
    cont = 1000
    mat1 = np.random.rand(ndim,ndim)
    matre = np.zeros((ndim,ndim))
    for i in range(cont):
        for j in range(cont):
            for k in range(cont):
                matre = np.dot(mat1,mat1)
    return matre

t1 = time.time()
# print(f1())   # 计算循环求和
print(f2())
t2 = time.time()
print('Timing cost is(s):  ',t2 - t1)
```

![png](/assets/images/Julia/p-mat1.png)

这里的测试表明这个循环加速过程对于矩阵的乘法也是同样适用的，服务器的所有核都被调用起来执行.

这里我想强调一下，刚开始我的矩阵相乘用的是`mat1*mat1`,我发现这种计算方式在利用jit进行并行加速的时候是失败的，如果采用numpy库的`np.dot(mat1,mat1)`结果就是正确的，可以很好的进行并行加速。
{:.warning}

# nodal-line 杂质计算程序
下面的程序是我想用来重复[Impurity-induced resonant states in topological nodal-line semimetals](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.205119)这篇文章的，程序的正确与否还正在检验中，不过确定的是在利用jit进行函数闭包之后，在执行速度上确实是由很明显的提升。
```python
from math import *   # 引入sqrt(), pi, exp等
import numpy as np
from numba import jit
import cmath 
import time
#=============================================
#@jit(nopython = True, parallel = True)
@jit(parallel = True)
def hamset(kx,ky,on):
    hn = on*2
    tx = 1.0
    ty = 1.0
    tz = 1.0
    lamz = 1.0
    ham = np.zeros((hn,hn)) * (1 + 0j)
    for i in range(40):
        ham[i,i + on] = tx*cos(kx) + ty*cos(ky)
        if(i != on - 1):
            ham[i,on + i + 1] = tz/2.0
        if(i != 0):
            ham[i,on + i - 1] = tz/2.0
        if(i != on - 1):
            ham[i,i + 1] = lamz/(2*1j)
        if(i != 0):
            ham[i,i - 1] = -lamz/(2*1j)
        ham[on + i,i] = tx*cos(kx) + ty*cos(ky)
        if(i != on - 1):
            ham[on + i,i + 1] = tz/2.0
        if(i != 0):
            ham[on + i,i - 1] = tz/2.0
        if(i != on - 1):
            ham[on + i,on + i + 1] = lamz/(2*1j)
        if(i != 0):
            ham[on + i,on + i - 1] = -lamz/(2*1j)
    return ham
#===============================================
#@jit(nopython = True, parallel = True)
@jit(parallel = True)
def spectra(kx,ky,on):
    gam = 0.01
    ham = hamset(kx,ky,on)
    eigval, eigvec = np.linalg.eig(ham)
    hn = len(ham)
    en = 200
    omglist = np.linspace(-cmath.pi, cmath.pi, en)
    G0 =  np.zeros((hn,hn,en)) * (1 + 0j)
    Gr0 =  np.zeros((hn,hn,en)) * (1 + 0j)
    G0r =  np.zeros((hn,hn,en)) * (1 + 0j)
    G0rr =  np.zeros((hn,hn,en)) * (1 + 0j)
    # 在所有格点上计算谱函数
    for m in range(hn):
        for n in range(hn):
            for kk in range(en):
                re1 = 0 + 0j
                for j in range(hn):
                    re1 += eigvec[n,j]*np.conj(eigvec[m,j])/(omglist[kk] - np.real(eigval[j]) + gam*1j)
                G0[m,n,kk] = re1
                Gr0[m,n,kk] = re1*cmath.exp(kx*1j)   
                G0r[m,n,kk] = re1*cmath.exp(-kx*1j)
                G0rr[m,n,kk] = re1*cmath.cos(ky)
    return G0,Gr0,G0r,G0rr        
#==============================================================
@jit(parallel = True)
def GreenFun(on,kn):
    kxlist = np.linspace(-cmath.pi, cmath.pi, kn)
    kylist = np.linspace(-cmath.pi, cmath.pi, kn)
    a1,a2,a3,a4 = spectra(0.1,0.1,on)
    l1,l2,l3 = np.shape(a1)
    G0 =  np.zeros((l1,l2,l3)) * (1 + 0j)
    Gr0 =  np.zeros((l1,l2,l3)) * (1 + 0j)
    G0r =  np.zeros((l1,l2,l3)) * (1 + 0j)
    G0rr =  np.zeros((l1,l2,l3)) * (1 + 0j)
    for ky in kylist:
        for kx in kxlist:
            a1,a2,a3,a4 = spectra(kx,ky,on)
            G0 += a1
            Gr0 += a2
            G0r += a3
            G0rr += a4
    # 对积分后的量乘以积分步长
    G0 = G0 * 2 * pi/kn
    Gr0 = Gr0 * 2 * pi/kn
    G0r = G0r * 2 * pi/kn
    G0rr = G0rr * 2 * pi/kn
    return G0,Gr0,G0r,G0rr
#==================================================================================
# @jit(parallel = True)
@jit
def Tmat(on,kn):
    v1 = 5;v2 = 10;v3 = 15;v4 = 20;v5 = 25;v6 = 30
    a1,a2,a3,a4 = spectra(0.1,0.1,on)
    d1,d2,d3 = np.shape(a1)  # d3得到的是omega撒点的个数,d1和d2得到是开边界格点数目*2
    # omglist = np.linspace(-cmath.pi, cmath.pi, d3) # omega的撒点取值
    one = np.identity(d1)
    U1 = np.zeros((d1,d2)) * (1.0 + 0j)
    U2 = np.zeros((d1,d2)) * (1.0 + 0j)
    U3 = np.zeros((d1,d2)) * (1.0 + 0j)
    U4 = np.zeros((d1,d2)) * (1.0 + 0j)
    U5 = np.zeros((d1,d2)) * (1.0 + 0j)
    U6 = np.zeros((d1,d2)) * (1.0 + 0j)
    ldos1 = np.zeros((d3,1)) * (1.0 + 0j)
    ldos2 = np.zeros((d3,1)) * (1.0 + 0j)
    ldos3 = np.zeros((d3,1)) * (1.0 + 0j)
    ldos4 = np.zeros((d3,1)) * (1.0 + 0j)
    ldos5 = np.zeros((d3,1)) * (1.0 + 0j)
    ldos6 = np.zeros((d3,1)) * (1.0 + 0j)
    U1[0,0] = v1
    U1[on,on] = v1
    U2[0,0] = v2
    U2[on,on] = v2
    U3[0,0] = v3
    U3[on,on] = v3
    U4[0,0] = v4
    U4[on,on] = v4
    U5[0,0] = v5
    U5[on,on] = v5
    U6[0,0] = v6
    U6[on,on] = v6
    G0,Gr0,G0r,G0rr = GreenFun(on,kn)
    for momg in range(d3):
        T1 = np.dot(np.linalg.inv(one - U1 * G0[:,:,momg]), U1)
        T2 = np.dot(np.linalg.inv(one - U2 * G0[:,:,momg]), U2)
        T3 = np.dot(np.linalg.inv(one - U3 * G0[:,:,momg]), U3)
        T4 = np.dot(np.linalg.inv(one - U4 * G0[:,:,momg]), U4)
        T5 = np.dot(np.linalg.inv(one - U5 * G0[:,:,momg]), U5)
        T6 = np.dot(np.linalg.inv(one - U6 * G0[:,:,momg]), U6)

        Grr1 = G0rr[:,:,momg] + np.dot(np.dot(Gr0[:,:,momg] ,T1) , G0r[:,:,momg])
        Grr2 = G0rr[:,:,momg] + np.dot(np.dot(Gr0[:,:,momg] ,T2) , G0r[:,:,momg])
        Grr3 = G0rr[:,:,momg] + np.dot(np.dot(Gr0[:,:,momg] ,T3) , G0r[:,:,momg])
        Grr4 = G0rr[:,:,momg] + np.dot(np.dot(Gr0[:,:,momg] ,T4) , G0r[:,:,momg])
        Grr5 = G0rr[:,:,momg] + np.dot(np.dot(Gr0[:,:,momg] ,T5) , G0r[:,:,momg])
        Grr6 = G0rr[:,:,momg] + np.dot(np.dot(Gr0[:,:,momg] ,T6) , G0r[:,:,momg])

        ldos1[momg] =  -(1/cmath.pi) * np.imag(Grr1[0,0] + Grr1[on,on])
        ldos2[momg] =  -(1/cmath.pi) * np.imag(Grr2[0,0] + Grr2[on,on])
        ldos3[momg] =  -(1/cmath.pi) * np.imag(Grr3[0,0] + Grr3[on,on])
        ldos4[momg] =  -(1/cmath.pi) * np.imag(Grr4[0,0] + Grr4[on,on])
        ldos5[momg] =  -(1/cmath.pi) * np.imag(Grr5[0,0] + Grr5[on,on])
        ldos6[momg] =  -(1/cmath.pi) * np.imag(Grr6[0,0] + Grr6[on,on])
    return ldos1,ldos2,ldos3,ldos4,ldos5,ldos6
#=================================================================================
@jit
def dataSave(on,kn):
    d1,d2,d3,d4,d5,d6 = Tmat(on,kn)
    en = len(d1)
    olist = np.linspace(-cmath.pi, cmath.pi, en).reshape((en,1))
    re = np.hstack((np.real(olist),np.real(d1),np.real(d2),np.real(d3),np.real(d4),np.real(d5),np.real(d6)))
    np.savetxt('ldos.dat',np.real(re),fmt='%.5f')

#==================================================================================
def main():
    on = 40  # open lattice size
    kn = 512  # k-point integration
    t1 = time.time()
    dataSave(on,kn)
    t2 = time.time()
    print('Timing cost is(s):  ',t2 - t1)
#========================================================
if __name__ == '__main__':
    main()

```