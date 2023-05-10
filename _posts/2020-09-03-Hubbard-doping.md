---
title: 电子掺杂，空穴掺杂与化学势的关系
tags: Fortran Superconductor
layout: article
license: true
toc: true
key: a20200831
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
在研究超导问题的时候，不可避免的总会遇到体系时电子掺杂还是空穴掺杂，因为这两种不同的掺杂所对应的性质是不同的，最近重新温习一下超导的相关知识，正好也利用一个具体的紧束缚近似模型来说明一下电子和空穴掺杂到底是怎么回事，以及如何的通过自洽方法，通过调整化学势来决定体系到底是电子掺杂还是空穴掺杂
{:.info}
<!--more-->
# 半满
## 实空间分析
在这里先采用李正中**<<固体理论>>**第十一章第九节的[Hubbard](https://en.wikipedia.org/wiki/Hubbard_model)模型来说明一下填充到底如何确定
$$
H=-t\sum_{<ij>,\sigma}(C_{i\sigma}^\dagger C_{j\sigma}+h.c)+U\sum_in_{i\uparrow}n_{i\downarrow}
$$
最后一项就是相同格点之内的不同自旋电子之间的相互作用，熟悉强关联的一定对这个模型非常熟悉了，而正好也就是这一项，使得本来可以利用紧束缚近似完整描述的模型，变成了强关联体系，从而引起了一系列有趣的现象。

在这里简单的令$U=0$，这样剩余的部分就是一个简单的紧束缚近似模型，学过固体物理就很清楚，如果你有N个格点的话，利用Pauli不相容原理，每个格点最多放两个电子，那么整个体系的话就最多可以放入2N个电子，这种情况就是**满填充**。那半满就很容易理解了，平均计算下来，每个格点上填充1个电子，那么整个体系存在N个电子，就是所谓的**半满**。也就是说格点i上不同自旋电子数的和为1，如下所示
$$
n=<n_{i\uparrow}>+<n_{i\downarrow}>=1
$$
以半满情况下的费米面为基本零点，就类似于经典力学种选水平面为重力势能零点一样的思想，也就是说半满情况的化学势**=0**，那么如果化学势往上调，则就是向体系种加入更多的电子，这个时候就对应着电子掺杂，每个格点上的平均粒子数会大于1。相反的如果把化学势向下调，那么就是从体系中拿走更多的电子，这时候每个格点上的平均粒子数就会小于1，这也就是空穴掺杂。对于掺杂这件事请，从理论的角度上看确实挺简单，但是在实验上你总是可以通过对某一种元素进行化学替代，从而实现电子或者空穴掺杂。
## 动量空间分析
熟悉紧束缚近似就一定明白，实空间和动量空间种的Hamiltonian也就是通过Fourier变换相互联系的。而在动量空间中你只不过是没有了格点的概念，不过它对应到动量空间中正好就是动量k，而且k点的数量和你变换到实空间中格点的数量也是一一对应的。

所以在动量空间中计算的时候就是对于每个动量格点k，计算其对应的能量大小。在上面也说过，以半满时候的费米面为界限，此时的化学势=0，那么如果某个动量点k其对应的能量大于0，那么它此时就是空的，因为费米面一下的状态才是被电子完全占据的，那么如果某个k点对应的能量是小于0的，那么就说明这个点是占据态，那么也就是说这个k点上对应的电子数为2，空态对应的电子数自然就是0，所以空态就不用去考虑了。

根据上面的分析，遍历所有的k点之后，你总可以确定哪些k点对应的是占据态，哪些点对应的是空态，将所有占据态电子数进行求和之后，再平均到你遍历时的每一个k点上，那么结果就是每个k点对应的电子占据数是多少。

>这里可能会有些不清楚，最后是将总的粒子数平均到所有的k点上，包括那些是空态的k点上，这样才是正确的。其实这里和实空间中在每个格点上计算平均粒子数是相同的逻辑

# 实例计算
## 化学势自洽
```fortran
subroutine chem()
use pub
integer m1,m2,m3
real kx,ky,kz,re1,re2,dev,filling
real,external::eng
do while(.true.)
    filling = 0
    do m1 = -kn,kn
        do m2 = -kn,kn
            do m3 = -kn,kn
                kx = pi*m1/kn
                ky = pi*m2/kn
                kz = pi*m3/kn
                re1 = eng(kx,ky,kz)
                if(re1 < 0)then
                    filling = filling + 2.0  ! 费米面以下满填充，每个态存在两个电子
                end if
            end do
        end do
    end do
    filling = filling/(2.0*kn + 1.0)**3
    !------------------------------
    !通过doping来控制掺杂浓度，和上面的保持一致，以半满填充为分解线，大于1就是电子掺杂，小于1就是空穴掺杂
    dev = filling - (1.0 - doping)  ! 检验计算得到的k点上的平剧粒子数与掺杂要求是否一致
    if(abs(dev) < err)then
        mu = mu
        ! call band(0.0)
        exit
    else
        if(dev>0)then
            mu = mu - dev*0.2   !  计算所得填充比要求高，则下调化学势
        else
            mu = mu + dev*0.2   !  计算所得填充比要求低，则上调化学势
        end if
    end if
    open(12,file="mu.dat",access="append")
    write(12,*)mu,filling,dev
    close(12)
end do
return
end subroutine chem
```
## 紧束缚能带
```fortran
real function eng(kx, ky, kz)
use pub
real kx,kz,ky
eng = -tp1*2.0*(cos(kx) + cos(ky)) - tp2*4.0*cos(kx)*cos(ky) - tp3*2.0*(cos(2.0*kx) + cos(2.0*ky))-tv1*cos(kz)*((cos(kx) - cos(ky))**2.0)/4.0 -tv2*cos(2.0*kz)*(cos(kx) - cos(ky))**2.0/4.0 - mu
end function eng
```

# 参考
上面紧束缚的能带是来源于这篇文章[Spin excitations in nickelate superconductors
](https://arxiv.org/abs/1910.05757v1)，我也利用上面的程序去重复了文章中提到的掺杂浓度和化学式之间的对应值，感兴趣可以自己试试。

关于实空间中格点上计算掺杂和化学势虽然逻辑和动量空间中的类似，但是程序上会有一些不同的地方，之后我会单独整理一个Blog来详细的再说明实空间的的这种计算要如何实现
{:.info}

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)