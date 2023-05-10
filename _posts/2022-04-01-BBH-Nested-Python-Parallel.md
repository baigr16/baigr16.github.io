---
title: BBH Nested Wilson loop计算(Python 并行)
tags: Python Topology 
layout: article
license: true
toc: true
key: a20220401
pageview: true
cover: /assets/images/python/Nested-1.png
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
这里就使用`Python`写了一下Nested Wilson loop的并行程序，因为前面考虑的体系其实就是一个$4\times 4$的哈密顿量，如果体系变大以及计算撒点数量变大的时候，`Python`还是会有点慢，这里就给一个并行的版本，并顺便看一下`Python`怎么实现进行并行。
{:.info}
<!--more-->
# Nested Wilson Loop
关于Nested Wilson loop的计算可以参考[Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators
](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)这篇文章，这里主要说一下自己在学习计算的时候踩过的坑。

与计算Wilson loop相同，这里最主要的仍然是找到一个Wannier band basis，也就是文章的中的公式

$$\rvert w_{x,\mathbf{k}}^j\rangle=\sum_{n=1}^\text{Nocc}\rvert u^n_\mathbf{k}\rangle[v_{x,\mathbf{k}}^j]^n$$

其实在做计算的时候，最让人困扰的不过是公式中的一大堆符号对应的到底是什么，这里就来讲这个公式拆解开，一步一步的讲清楚里面的含义。这里假设你已经知道为什么要计算Nested Wilson loop，
我在这里就简单的阐述一下。首先要是体系的Wilson loop计算对应的Wannier哈密顿量的能带是有能隙的，也就是说你的体系是4带模型，那么当占据态是2条能带的时候，每个占据态能带会对应着一个Wannier center，
比如BHZ模型的两条占据态能带对应的Wannier band就是相互交叉的，而且因为Wilson loop与边界态之间的拓扑等价性，TI是有边界态的，所以其对应的Wilson loop在形状上就与边界态类似。而对于高阶拓扑相，
首先就是要使得边界态打开能隙，那么相对应的就是其Wilson loop计算得到的Wannier center随着某个动量参数的演化是不会相互交叉的，这一点在上面BBH模型中已经计算过了，所以此时就可以对某一个单独的Wannier band
计算它的Nested Wilson loop，所以首先第一步就是必须要明白什么样的情况下，是需要计算体系的Nested Wilson loop。

这里的$\sum_{n=1}^\text{Nocc}$不用讲太多，是需要对占据态进行求和，但是这个$n$其实表示的只是说哈密顿量的占据态，也就是说对于$\rvert u^n_\mathbf{k}\rangle$而言，这是哈密顿量的占据态波函数，$n$表示占据态其实是对$\rvert u^n_\mathbf{k}\rangle$
而言的，虽然$[v_{x,\mathbf{k}}^j]^n$中同样存在这个$n$，但是在这个地方$n$代表的不是占据态，在这里$j$代表的才是你选择的是哪一个Wannier band来计算其对应的Nested Wilson loop，也就是这里$j$代表的是你选择的是那个占据的Wannier band，而$n$在这里
表示的是一个Wannier哈密顿量本征矢量的第$n$个分量。假如$H_\text{Wann}$是Wannier哈密顿量，其满足


$$H_\text{Wann}\rvert v_\mathbf{k}\rangle=E\rvert v_\mathbf{k}\rangle$$

那么这里的$[v_{x,\mathbf{k}}^j]^n$表示的就是这个本征矢量的第$n$个分量，$j$则表示你选择是哪个本征值对应的本征矢量，也就是选择了哪一个Wannier band。这里的$x$则表示你在做Wilson loop的时候，是沿着哪个方向进行的，即就是讲上面公式中的$H_\text{Wann}$替换成你
构建的那个Wilson loop的哈密顿量就可以。

至于$\rvert u^n_\mathbf{k}\rangle$就很简单了，它表示的就是你的哈密顿量的本征态，当然了在计算的时候，还是要选择正确的占据态才可以。下面直接上代码，在其中同样做了注释
```python
import numpy as np
# import matplotlib.pyplot as plt
import os,time
import multiprocessing as mp # 进程并行
# from functools import partial # 固定参数
def Pauli():
    s0 = np.array([[1,0],[0,1]])
    sx = np.array([[0,1],[1,0]])
    sy = np.array([[0,-1j],[1j,0]])
    sz = np.array([[1,0],[0,-1]])
    return s0,sx,sy,sz  
#-------------------------------------------------------------
def hamset(kx,ky,m0):
    hn = 4
    gamx = 0.5  
    lamx = 1.0  
    gamy = gamx
    lamy = lamx
    xsyb1 = 0.000000000000    
    xsyb2 = 1.1
    ysyb1 = 0.000000000000    
    ysyb2 = 1.0
    ham = np.zeros((hn, hn),dtype = complex)
    ham[0, 0] = xsyb1
    ham[1, 1] = ysyb1
    ham[2, 2] = ysyb1
    ham[3, 3] = xsyb1
    ham[0, 2] = (gamx + lamx * np.exp(1j * kx)) * ysyb2
    ham[1, 3] = gamx + lamx * np.exp(-1j * kx)
    ham[0, 3] = gamy + lamy * np.exp(1j * ky)
    ham[1, 2] = (-gamy - lamy * np.exp(-1j * ky)) * xsyb2
    ham[2, 0] = np.conj(ham[0, 2])
    ham[3, 1] = np.conj(ham[1, 3])
    ham[3, 0] = np.conj(ham[0, 3])
    ham[2, 1] = np.conj(ham[1, 2])
    return ham
#-------------------------------------------------------------
def Wilson_kx(kx,m0,nk):
    hn = hamset(0,0,m0).shape[0] # 获取哈密顿量的维度，方便构建占据态
    Nocc = int(hn/2)
    wave = np.zeros((hn,hn,nk),dtype = complex) # 存储所有ki点上的本征波函数
    wanvec = np.zeros((Nocc,Nocc),dtype = complex)
    klist = np.linspace(-np.pi, np.pi, nk)
    Wan = np.eye(Nocc,dtype = complex)
    F = np.zeros((Nocc, Nocc), dtype = complex) # Wannier Hamiltonian
    # 构建沿着y方向的Wilson loop，此时给定一个kx值进行一次Wilson loop计算
    ix = 0
    for ky in klist: # 计算沿着kx方向的Wilson loop
        eigval, eigvec = np.linalg.eigh(hamset(kx, ky, m0)) # 求解哈密顿量的本征矢量
        wave[:,:,ix] = eigvec[:,:] # 存储所有的本征波函数,用来后面计算Wilson loop
        ix += 1
    wave[:,:,nk - 1] = wave[:,:,0] # 首尾波函数相同,消除规范
    # 利用沿着ky方向计算的波函数来构建Wannier hamiltonian
    for i0 in range(nk - 1):
        for i1 in range(Nocc): # 直接通过循环,只遍历占据态的波函数，在设定化学势为零的时候，占据态是能带数的一半
            for i2 in range(Nocc): # 不同占据态在相邻kx格点上的波函数交叠
                F[i1,i2] = np.dot(wave[:,i1,i0 + 1].transpose().conj(),wave[:,i2,i0])
        U,s,V = np.linalg.svd(F)
        F = np.dot(U,V)
        Wan = np.dot(F, Wan) 
    eigval, eigvec = np.linalg.eig(Wan) # 这里Wannier 哈密顿量并不是厄米的，所以求解得到的本征值的顺序并不一定就是按顺序排列的
    mux = np.log(eigval)/(2*np.pi*1j) # 从Wannier哈密顿量中计算Wannier center
    for i0 in range(Nocc):
        wanvec[:,i0] = eigvec[:, np.argsort(np.real(mux))[i0]]
    return wanvec   
#---------------------------------------------------------------------------------------------------------
def Wilson_ky(ky,m0,nk):
    hn = hamset(0,0,m0).shape[0] # 获取哈密顿量的维度，方便构建占据态
    Nocc = int(hn/2)
    wave = np.zeros((hn,hn,nk),dtype = complex) # 存储所有ki点上的本征波函数
    klist = np.linspace(-np.pi, np.pi, nk)
    wanvec = np.zeros((Nocc,Nocc),dtype = complex)
    # 构建沿着y方向的Wilson loop，此时给定一个kx值进行一次Wilson loop计算
    ix = 0
    for kx in klist: # 计算沿着kx方向的Wilson loop
        eigval, eigvec = np.linalg.eigh(hamset(kx, ky, m0)) # 求解哈密顿量的本征矢量
        wave[:,:,ix] = eigvec[:,:] # 存储所有的本征波函数,用来后面计算Wilson loop
        ix += 1
            # 首尾波函数相同,消除规范
    wave[:,:,nk - 1] = wave[:,:,0] 
    # 利用沿着ky方向计算的波函数来构建Wannier hamiltonian
    Wan = np.eye(Nocc,dtype = complex)
    F = np.zeros((Nocc, Nocc), dtype = complex) # Wannier Hamiltonian
    for i0 in range(nk - 1):
        for i1 in range(Nocc): # 直接通过循环,只遍历占据态的波函数
            for i2 in range(Nocc):
                F[i1,i2] = np.dot(wave[:,i1,i0 + 1].transpose().conj(),wave[:,i2,i0])
        U,s,V = np.linalg.svd(F)
        F = np.dot(U,V)
        Wan = np.dot(F, Wan) 
    eigval, eigvec = np.linalg.eig(Wan) # 这里Wannier 哈密顿量并不是厄米的，所以求解得到的本征值的顺序并不一定就是按顺序排列的
    mux = np.log(eigval)/(2*np.pi*1j) # 从Wannier哈密顿量中计算Wannier center
    for i0 in range(Nocc):
        wanvec[:,i0] = eigvec[:, np.argsort(np.real(mux))[i0]]
    return wanvec 
#---------------------------------------------------------------------------------------------------------
def Nested_Wilson_loop_kx(m0,nk):
    # nk = 10 # 计算Nested Wilson loop时撒点的数目
    hn = hamset(0,0,m0).shape[0] # 获取哈密顿量的维度，方便构建占据态
    Nocc = int(hn/2) # 占据态能带数目
    klist = np.linspace(-np.pi, np.pi, nk)
    wave = np.zeros((hn,hn,nk),dtype = complex)
    wann_vec = np.zeros((nk,Nocc,Nocc),dtype = complex)
    pmulist = np.zeros((nk,int(Nocc/2)),dtype = float)

    for iy in range(nk):
        ky = klist[iy]
        wann_vec[iy,:,:] = Wilson_ky(ky,m0,nk) # 在固定ky的情况下，计算沿着kx方向的Wilson loop并得到对应的本征矢量

    for ix in range(nk):
        kx = klist[ix]
        for iy in range(nk): # 沿着一个方向遍历动量
            ky = klist[iy]
            eigval,eigvec = np.linalg.eigh(hamset(kx,ky,m0)) # 求解哈密顿量的本征矢量和本征值
            if ky != np.pi:
                wave[:,:,iy] = eigvec[:,:] # 将沿着ky方向的所有本征波函数存储
            else:
                wave[:,:,nk - 1] = wave[:,:,0] # 在边界上保证波函数首尾相接
        #------------------------------------------------------------------
        wmu =np. zeros((int(Nocc/2),hn,nk),dtype = complex) 
        for i3 in range(int(Nocc/2)):
            for iy in range(nk):
                for i4 in range(Nocc):
                        wmu[i3,:,iy] += wave[:,i4,iy] * wann_vec[iy,i4,i3]
                wmu[i3,:,nk - 1] = wmu[i3,:,0] # 首尾相接
        #--------------------------------------
        wan = np.zeros((int(Nocc/2),int(Nocc/2)),dtype = complex)  # 对确定的Wannier sector中的每一个Wannier 能带计算Wilson loop
        F0 = np.zeros((int(Nocc/2),int(Nocc/2)),dtype = complex)
        for i4 in range(int(Nocc/2)):
            wan[i4,i4] = 1
        for iy in range(nk - 1):
            for i3 in range(int(Nocc/2)):
                for i2 in range(int(Nocc/2)):
                    F0[i3,i2] = np.dot(wmu[i3,:,iy + 1].transpose().conj(),wmu[i2,:,iy]) # 在新的Wannier basis下面构建Wilson loop，也就是计算Nested Wilson loop
            U,s,V = np.linalg.svd(F0)
            F0 = np.dot(U,V)
            wan = np.dot(F0,wan)
        val,vec = np.linalg.eig(wan)
        for i0 in range(len(val)):
            val[i0] = np.angle(val[i0])/(2 * np.pi)
            # if val[i0] < 1:
            #     val[i0] += 1
        pmulist[ix,:] = np.real(val)
    return pmulist
#---------------------------------------------------------------------------------------------------------
def Nested_Wilson_loop_ky(m0,nk):
    # nk = 10 # 计算Nested Wilson loop时撒点的数目
    hn = hamset(0,0,m0).shape[0] # 获取哈密顿量的维度，方便构建占据态
    Nocc = int(hn/2) # 占据态能带数目
    klist = np.linspace(-np.pi, np.pi, nk)
    wave = np.zeros((hn,hn,nk),dtype = complex)
    wann_vec = np.zeros((nk,Nocc,Nocc),dtype = complex)
    pmulist = np.zeros((nk,int(Nocc/2)),dtype = float)

    for ix in range(nk):
        kx = klist[ix]
        wann_vec[ix,:,:] = Wilson_kx(kx,m0,nk) # 在固定ky的情况下，计算沿着kx方向的Wilson loop并得到对应的本征矢量

    for iy in range(nk):
        ky = klist[iy]
        for ix in range(nk): # 沿着一个方向遍历动量
            kx = klist[ix]
            eigval,eigvec = np.linalg.eigh(hamset(kx,ky,m0)) # 求解哈密顿量的本征矢量和本征值
            if ky != np.pi:
                wave[:,:,ix] = eigvec[:,:] # 将沿着ky方向的所有本征波函数存储
            else:
                wave[:,:,nk - 1] = wave[:,:,0] # 在边界上保证波函数首尾相接
            #------------------------------------------------------------------
        wmu =np. zeros((int(Nocc/2),hn,nk),dtype = complex) 
        for i3 in range(int(Nocc/2)):
            for ix in range(nk):
                for i4 in range(Nocc):
                        wmu[i3,:,ix] += wave[:,i4,iy] * wann_vec[iy,i4,i3]
                wmu[i3,:,nk - 1] = wmu[i3,:,0] # 首尾相接
        #--------------------------------------
        wan = np.zeros((int(Nocc/2),int(Nocc/2)),dtype = complex)  # 对确定的Wannier sector中的每一个Wannier 能带计算Wilson loop
        F0 = np.zeros((int(Nocc/2),int(Nocc/2)),dtype = complex)
        for i4 in range(int(Nocc/2)):
            wan[i4,i4] = 1
        for ix in range(nk - 1):
            for i3 in range(int(Nocc/2)):
                for i2 in range(int(Nocc/2)):
                    F0[i3,i2] = np.dot(wmu[i3,:,ix + 1].transpose().conj(),wmu[i2,:,ix]) # 在新的Wannier basis下面构建Wilson loop，也就是计算Nested Wilson loop
            U,s,V = np.linalg.svd(F0)
            F0 = np.dot(U,V)
            wan = np.dot(F0,wan)
        val,vec = np.linalg.eig(wan)
        for i0 in range(len(val)):
            val[i0] = np.angle(val[i0])/(2 * np.pi)
            # if val[i0] < 1:
            #     val[i0] += 1
        pmulist[iy,:] = np.real(val)
    return pmulist
#-------------------------------------------------------
def Nested1(m0):
    cont = 1
    nk = 100 # 控制计算撒点数目
    klist = np.linspace(-np.pi, np.pi, nk)
    x1 = Nested_Wilson_loop_kx(m0,nk)
    x2 = Nested_Wilson_loop_ky(m0,nk)
    re1 = 0
    re2 = 0
    for i0 in range(nk):
        re1 += np.mod(np.sum(x1[i0,:]),1)
        re2 += np.mod(np.sum(x2[i0,:]),1)
    re1 = re1/nk
    re2 = re2/nk
    # f1 = "Nested-kx-" + str(cont) + ".dat"
    # np.savetxt(f1,x1,fmt="%11.6f")
    # f2 = "Nested-ky-" + str(cont) + ".dat"
    # np.savetxt(f2,x2,fmt="%11.6f")
    return m0,re1,re2,2*re1*re2
#------------------------------------------------------------------
def main():
    num = 40 # 控制质量项变化时的撒点数量
    m0list = np.linspace(0,3, num)
    print("Number of process :%d"%(mp.cpu_count()))
    pool = mp.Pool(mp.cpu_count())
    #pool = mp.Pool(3)
    result = pool.map(Nested1,m0list)  # 同步提交
    re2 = np.mat(list(map(list,result)))
    np.savetxt("te-polar-100.dat",re2,fmt="%11.6f")
#------------------------------------------------------------------
if __name__ == '__main__':
    os.chdir(os.getcwd()) # 更改工作目录到当前文件夹
    tstart = time.time() # 获取系统时间
    main()
    tend = time.time()
    print("Consted time is %.5f" % (tend - tstart))
```
# Python 进程并行
```python
def main():
    num = 40 # 控制质量项变化时的撒点数量
    m0list = np.linspace(0,3, num)
    print("Number of process :%d"%(mp.cpu_count()))
    pool = mp.Pool(mp.cpu_count())
    #pool = mp.Pool(3)
    result = pool.map(Nested1,m0list)  # 同步提交
    re2 = np.mat(list(map(list,result)))
    np.savetxt("te-polar-100.dat",re2,fmt="%11.6f")
```
其实`Python`改并行还是比较方便的，只要将一个循环的变量存储为`list`，然后`map`到不同的进程上就可以
```python
m0list = np.linspace(0,3, num)
pool = mp.Pool(mp.cpu_count()) # cpu有所少进程，就开启多少进程，当然在这里你也可以指定进程数量
result = pool.map(Nested1,m0list)  # 这一步就是将一个循环过程分发到不同的进程上，最后所有的结果都会返回出来
```
上面的程序是比较通用的，要想计算别的模型，就只需要将哈密顿量替换一下即可。

# 绘图程序
这里就用`python`画图了，给一个绘图脚本
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
config = {
"font.size": 30,
"mathtext.fontset":'stix',
"font.serif": ['SimSun'],
}
rcParams.update(config) # Latex 字体设置
#---------------------------------------------------------
def scatterplot1(cont):
    da1 = "Nested-kx-" + str(cont) + ".dat"
    da2 = "Nested-ky-" + str(cont) + ".dat"
    picname = "Nested-" + str(cont) + ".png"
    os.chdir(os.getcwd())# 确定用户执行路径
    x0 = []
    y0 = []
    with open(da1) as file:
        da = file.readlines()
        for f1 in da:
            if len(f1) > 3:
                ldos = [float(x) for x in f1.strip().split()]
                x0.append(ldos[0])
                y0.append(ldos[1])
    y0 = np.array(y0)
    plt.scatter(x0, y0, s = 20, color = 'lightskyblue', label = "$p_y^{v_x^\pm}(k_x)$")
    x1 = []
    y1 = []
    with open(da2) as file:
        da = file.readlines()
        for f1 in da:
            if len(f1) > 3:
                ldos = [float(x) for x in f1.strip().split()]
                x1.append(ldos[0])
                y1.append(ldos[1])
    y1 = np.array(y1)
    # print(y0)
    # sc = plt.scatter(x0, y0, c = z1, s = 2,vmin = 0, vmax = 1, cmap="magma")
    plt.scatter(x1, y1, s = 20, color = 'deeppink', label = "$p_x^{v_y^\pm}(k_y)$")
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 25,
             }
    plt.xlim(0,1)
    plt.ylim(-1,1)
    plt.xlabel("$k_x(k_y)/\pi$",font2)
    # plt.ylabel("",font2)
    plt.yticks([-1,-0.5,0.,0.5,1],fontproperties='Times New Roman', size = 25)
    plt.xticks([-1,0,1],fontproperties='Times New Roman', size = 25)
    plt.legend(loc = 'upper left', bbox_to_anchor=(0.5,0.5), shadow = True, prop = font2, markerscale = 4)
    # plt.text(x = 0.6,y = 0.7,s = 'MCM', fontdict=dict(fontsize=20, color='black',family='Times New Roman'))
    # plt.text(x = 0.1,y = 0.7,s = 'NSC', fontdict=dict(fontsize=20, color='black',family='Times New Roman'))
    # plt.vlines(x = 0.4, ymin = -1, ymax = 1,lw = 3.0, colors = 'black', linestyles = '--')
    plt.savefig(picname, dpi = 600, bbox_inches = 'tight')
#---------------------------------------------------------
def main():
    for i0 in range(1,2):
        scatterplot1(i0) 
#---------------------------------------------------------
if __name__=="__main__":
    main()
```

![png](/assets/images/python/Nested-1.png)

# 参考
- 1.[Electric multipole moments, topological multipole moment pumping, and chiral hinge states in crystalline insulators
](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.245115)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)