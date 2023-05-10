---
title: 矩阵对角化与Bogoliubov对角化的联系
tags: Study Method
layout: article
license: true
toc: true
key: a20200925b
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
在[前面](https://yxli8023.github.io/2020/09/25/MF-BdG.html)已经做过了对二次型哈密顿量进行Bogoliubov对角化的操作,但平时在做计算的时候,通常都是通过数值的方法来对哈密顿量矩阵进行对角化.在这里正好就将这两种方式求解之间的关系进行一下梳理.
{:.info}
<!--more-->
# 矩阵对角化
设M是一个任意的矩阵$(n\times n)$,假设它的本征值为$\lambda_1,\lambda_2,\dots,\lambda_n$,本征矢量为$\boldsymbol{V_1},\boldsymbol{V_2},\dots,\boldsymbol{V_n}$,这些本征矢量之间是相互正交的,形成线性无关的集合,利用这些本征矢量来构建一个矩阵

$$\boldsymbol{A} = [\boldsymbol{V_1},\boldsymbol{V_2},\dots,\boldsymbol{V_n}]$$

可以通过作用$\boldsymbol{A}$矩阵到$\boldsymbol{M}$上对其进行对角化

$$\boldsymbol{A}^{-1} \boldsymbol{M} \boldsymbol{A}=\left(\begin{array}{cccc}
\lambda_{1} & 0 & \cdots & 0 \\
0 & \lambda_{2} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & \lambda_{n}
\end{array}\right)$$

上面这也是一个相似变换过程,将这个过程反过来可以得到$\boldsymbol{M} = \boldsymbol{A}\boldsymbol{D}\boldsymbol{A^{-1}}$,这里$\boldsymbol{D}$是只有对角线元素的矩阵,每个对角元素正好对应着$\boldsymbol{M}$的本征值.

# 哈密顿量矩阵

$$\bar{H} \approx \sum_{k} \varepsilon_{k}\left(C_{k}^{+} C_{k}+C_{-k}^{+} C_{-k}\right)-\Delta \sum_{k}\left(C_{k}^{+} C_{-k}^{+}+C_{-k} C_{k}\right)+\Delta^{2} / V$$

对于这个哈密顿量,选取$\psi^\dagger = (C^\dagger_k,C_k,C^\dagger_{-k},C_{-k})$为基矢构建上述哈密顿量的矩阵形式,这个基矢的构建也就是在Nambu表象上进行的,也可以参考[这里](https://www.manep.ch/saasfee15/pdf/Black-Schaffer-1-2x2.pdf).

$$H = \psi^\dagger\mathcal{H}\psi=(C^\dagger_k,C_k,C^\dagger_{-k},C_{-k})\mathcal{H}_{4\times 4}\left(\begin{array}{l}
C_k&\\
C^\dagger_k&\\
C_{-k}&\\
C^\dagger_{-k}&
\end{array}\right)\label{ham}$$

$$\mathcal{H}=\left[\begin{array}{cccc}
\epsilon_k&0&0&-\Delta&\\
0&-\epsilon_k&\Delta&0&\\
0&\Delta&\epsilon_{-k}&0&\\
-\Delta&0&0&-\epsilon_{-k}&\\
\end{array}\right]$$

在研究超导问题的时候,基矢要比平时扩大一倍,哈密顿量是具有[粒子空穴对称性](https://physics.stackexchange.com/questions/86293/what-is-the-definition-of-particle-hole-symmetry-in-condensed-matter-physics),所以这里的空间是有冗余的.

已经有了矩阵之后就可以对其进行对角化了,当然这个过程是通过程序进行了,之后可以得到其对应的本征值和本征矢量,下面是一个简单的示例程序,在调用cheevd函数之后,Ham中存储的就是本征矢量,而w变量中存储的就是本征值,这里数值对角化的过程就完成了.
```fortran
! Author:YuXuanLi
! E-Mail:yxli406@gmail.com
! Article:Majorana Corner Modes in a High-Temperature Platform
! Doi:10.1103/PhysRevLett.121.096803
!==========================================
    module param
    implicit none
    integer xn,yn,nkx,nky,ne,nkxy
    parameter(nkx=50,nky=50,ne=1000)
    integer,parameter::N = 4
    complex,parameter::im = (0.,1.) !Imagine unit
    real,parameter::pi = 3.14159265358979
    complex Ham(N,N) ! Hamiltonian Matrix
    real mu ! Chemical Potential
    real tx,ty  ! hopping term energy
    real ax,ay  ! copule energy
    real m0  !Driac mass
    real kx,ky
    ! LAPACK PACKAGE PARAM
    integer::lda = N
    integer,parameter::lwmax = 2*N+N**2
    real,allocatable::w(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   ! at least 2*N+N**2
    integer lrwork    ! at least 1 + 5*N +2*N**2
    integer liwork   ! at least 3 +5*N
    integer info
    end module param
!========== PROGRAM START ==========================
    program sol
    use param
    integer m,l ! loop variales
    !================ Physics memory allocate =================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1+5*N+2*N**2))
    allocate(iwork(3+5*N))
    ! paramater set value
    m0 = 1.5   ! effective mass
    tx = 1.0   ! hopping energy in x direction
    ty = 1.0   ! hopping energy in y direction
    mu = 0.2   ! chemical potential
    ax = 1.0  ! couple energy in x direction
	  ay = 1.0  ! couple energy in y direction
    call band()
    stop
    end program
!===========================================================
    subroutine matrix_set()
    ! 构造哈密顿量矩阵
    use param
    Ham = 0
    !----------------------------------------------
    Ham(1,1) = m0 - tx*cos(kx) - ty*cos(ky) + mu
    Ham(2,2) = -(m0 - tx*cos(kx) - ty*cos(ky)) + mu
    Ham(3,3) = m0 - tx*cos(kx) - ty*cos(ky) + mu
    Ham(4,4) = -(m0 - tx*cos(kx) - ty*cos(ky)) + mu
    !-----------------------------------------------
    Ham(1,2) = ax*sin(kx) - im*ay*sin(ky) 
    Ham(2,1) = ax*sin(kx) + im*ay*sin(ky) 
    Ham(3,4) = -ax*sin(kx) - im*ay*sin(ky) 
    Ham(4,3) = -ax*sin(kx) + im*ay*sin(ky) 
    !--------------------------------------------------
    end subroutine matrix_set
!============================================================
    subroutine band()
    ! Evaluate the density of state in (x,y) position
    use param
    integer m,l,k,i! circle variable
    open(12,file="tiband.dat")
    !   (0,0)------>(pi,0)
    do k = 0,nky
        kx = pi*k/nkx  ! discrete wavevector in x direction
        !kx = 0
        ky = 0  ! discrete wavevector in y direction
        ! 不同的kx和ky重新进行哈密顿量矩阵的构造和求解
        call matrix_set() ! 新的kx和ky下重新填充矩阵，并求解对应本征值
        call eigSol()
        write(12,"(6f9.5)")kx,ky,(w(i),i=1,N)
    end do
    close(12)
    end subroutine band
!================= Hermitain Matrices solve ==============
      subroutine eigSol()
      use param
      integer m
      lwork = -1
      liwork = -1
      lrwork = -1
      call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
      lwork = min(2*N+N**2, int( work( 1 ) ) )
      lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
      liwork = min(3+5*N, iwork( 1 ) )
      call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
      if( info .GT. 0 ) then
            open(101,file="mes.txt",status="unknown")
            write(101,*)'The algorithm failed to compute eigenvalues.'
            close(101)
      end if
      open(100,file="eigval.dat",status="unknown")
      do m = 1,N
            write(100,*)m,w(m)
      end do
            close(100)
      return
      end subroutine eigSol

```

# Bogoliubov变换
下面还是对同一个哈密顿量,利用Bogoliubov变换对其进行对角化

$$\bar{H} \approx \sum_{k} \varepsilon_{k}\left(C_{k}^{+} C_{k}+C_{-k}^{+} C_{-k}\right)-\Delta \sum_{k}\left(C_{k}^{+} C_{-k}^{+}+C_{-k} C_{k}\right)+\Delta^{2} / V \label{bcs1}$$

引入新的粒子算符,关于这个方法的详细内容可以参考[这里](https://yxli8023.github.io/2020/09/25/MF-BdG.html).

$$\begin{aligned}
\alpha_{k} &=u_{k} C_{k}-v_{k} C_{-k}^{+}, \quad \alpha_{k}^{+}=u_{k} C_{k}^{+}-v_{k} C_{-k} \\
\alpha-k &=u_{k} C_{-k}+v_{k} C_{k}^{+}, \quad \alpha_{-k}^{+}=u_{k} C_{-k}^{+}+v_{k} C_{k}
\end{aligned}$$

对于这个新的算符,它同样满足费米子算符的对易关系

$$\begin{array}{c}
{\left[\alpha_{k}, \alpha_{k^{\prime}}^{+}\right]_{+}=u_{k} u_{k^{\prime}}\left[C_{k}, C_{k^{\prime}}^{+}\right]_{+}+v_{k} v_{k^{\prime}}\left[C_{-k}^{+}, C_{-k^{\prime}}\right]_{+}} \\
=\delta_{k k^{\prime}}\left(u_{k}^{2}+v_{k}^{2}\right)=\delta_{k k^{\prime}} \\
{\left[\alpha_{k}, \alpha_{-k^{\prime}}\right]_{+}=u_{k} v_{k^{\prime}}\left[C_{k}, C_{k^{\prime}}^{+}\right]_{+}-v_{k} u_{k^{\prime}}\left[C_{-k}^{+}, C_{-k^{\prime}}\right]_{+}} \\
=\delta_{k k^{\prime}}\left(u_{k} v_{k}-u_{k} v_{k}\right)=0 \\
{\left[\alpha_{-k}, \alpha_{-k^{\prime}}^{+}\right]_{+}=\delta_{k k^{\prime}}\left(u_{k}^{2}+v_{k}^{2}\right)=\delta_{k k^{\prime}}} \\
{\left[\alpha_{k}^{+}, \alpha_{-k^{\prime}}^{+}\right]_{+}=0}
\end{array}$$

根据上面准粒子算符$\alpha$的对易关系可以得到变换系数之间的一个关系

$$u_{k}^{2}+v_{k}^{2}=1\label{eq1}$$

因为变换是有两个系数的,此时只得到一个方程,还需要寻找另外一个方程,接下来将算符$C$利用准粒子算符表示出来

$$\left.\begin{array}{l}
C_{k}=u_{k} \alpha_{k}+v_{k} \alpha_{-k}^{+}, \quad C_{k}^{+}=u_{k} \alpha_{k}^{+}+v_{k} \alpha-k \\
C_{-k}=u_{k} \alpha_{-k}-v_{k} \alpha_{k}^{+}, \quad C_{-k}^{+}=u_{k} \alpha_{-k}^{+}-v_{k} \alpha_{k}
\end{array}\right\}$$

将这个算符的表达式回代到(\ref{bcs1})中,可以将其利用准粒子算符$\alpha$来表示

$$\begin{aligned}
\bar{H}=\sum_{k}\left\{\left[\varepsilon_{k}\left(u_{k}^{2}-v_{k}^{2}\right)+2 \Delta u_{k} v_{k}\right]\left(\alpha_{k}^{+} \alpha_{k}+\alpha_{-k}^{+} \alpha_{-k}\right)+\right.& \\
\left.\left[2 \varepsilon_{k} u_{k} v_{k}-\Delta\left(u_{k}^{2}-v_{k}^{2}\right)\right]\left(\alpha_{k}^{+} \alpha_{-k}^{+}+\alpha_{-k} \alpha_{k}\right)\right\}+& \\
\sum_{k}\left[2 \varepsilon_{k} v_{k}^{2}-2 u_{k} v_{k} \Delta\right]+\frac{\Delta^{2}}{V}
\end{aligned}\label{quasi}$$

Bogoliubov变换的主要目的就是要把哈密顿量变换成$\alpha^\dagger\alpha$这种形式,在上式中可以看到还存在$\alpha^\dagger\alpha^\dagger$这样的项,那么为了消除它们,令它们前面的系数为0.
{:.success}

所以在这里就可以得到关于$u_k,v_k$的另外一个方程

$$\Delta\left(u_{k}^{2}-v_{k}^{2}\right)=2 \varepsilon_{k} u_{k} v_{k}\label{eq2}$$

结合方程(\ref{eq1})和(\ref{eq2})可以求得变换的系数

$$u_{k}^{2}=\frac{1}{2}\left(1+\frac{\varepsilon_{k}}{\xi_{k}}\right), \quad v_{k}^{2}=\frac{1}{2}\left(1-\frac{\varepsilon_{k}}{\xi_{k}}\right)$$

这里的$\xi_k=\sqrt{\epsilon_k^2+\Delta^2}$,代表的就是超导态准粒子的激发能.

通过上面的一系列操作之后,哈密顿量变成了关于准粒子算符$\alpha$的对角形式

$$H\sim \xi_k\alpha_k^\dagger\alpha\label{eq5}$$

这里$\xi_k$也就是这个准粒子的本征能量,与前面数值对角化相比,这也就对应着哈密顿量的本征值.

# 相互联系

在前面首先知道矩阵$\boldsymbol{M} = \boldsymbol{A}\boldsymbol{D}\boldsymbol{A^{-1}}$,现在这个$\boldsymbol{M}$就是哈密顿量矩阵$\mathcal{H}$,那么同样的它的对角化为

$$\mathcal{H}=SDS^{-1}\label{eq3}$$

这里的S就是通过对角化得到的本征矢量构成的矩阵,也就是前面的$\boldsymbol{A}$,将(\ref{eq3})代入(\ref{ham})中,可以得到

$$H=\psi^\dagger SDS^{-1}\psi=(\psi S^\dagger)^\dagger D (S^\dagger\psi)\label{eq4}$$

这里D是对角矩阵,也就是本征值,将(\ref{eq4})和(\ref{eq5})进行对比就可以发现下列关系

$$
\alpha = S^\dagger\psi\\
\xi_k = D
$$

从这里饿哦们就可以很明确的看出矩阵对角化与Bogoliubove变换之间的联系了,矩阵对角化后的本征值就对应着Bogoliubov变换之后准粒子对角二次型前面的系数$\xi_k$.

这里因为$\mathcal{H}$在构建的时候,里面的数同样都是和k相关的参数,所以可以认为是不同的k对应着不同的$\mathcal{H}$,同样也就对应着不同的$D$,则可以将矩阵$D$在形式上也写作$D_k$,这样的话它个$\xi_k$之间的对应关系就很明朗了.
{:.warning}

至于算符之间的关系就更加明确了,准粒子算符$\alpha$和原始的费米子算符$C$之间通过一个幺正矩阵(厄米矩阵本征矢构成的矩阵是个幺正矩阵)联系,这也就和Bogoliubov变换时,准粒子算符由原始费米子算符通过系数组合联系起来,而且在执行这个准粒子算符构建过程的时候,这个变换本来就是幺正的,所有的内容到这里就变得完全自洽,而Bogoliubov变换的系数就可以通过矩阵对角化后得到的本征矢量矩阵$S$得到.
{:.success}

# 参考

- 1.固体理论(李正中)
- 2.[Majorana quasiparticles in condensed matter](https://arxiv.org/pdf/1711.00011.pdf)
- 3.[Topological Superconductivity ](https://www.manep.ch/saasfee15/pdf/Black-Schaffer-1-2x2.pdf)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)