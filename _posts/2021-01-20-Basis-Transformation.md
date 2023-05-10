---
title: BdG哈密顿量基矢修改
tags:  Superconductor
layout: article
license: true
toc: true
key: a20210120c
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
前面有两篇博客主要关注的不同的基矢下如何构建BdG形式的哈密顿量,既然只不过是基矢的变化,那么两个哈密顿量肯定就是等价的,这里就给一个具体的实例,来看看怎么把两个不同基矢下的哈密顿量进行相互的改写,并利用程序验证其正确性.
{:.info}
<!--more-->
# 模型
这里选用**参考(5)**文献中的哈密顿量

$$\begin{equation}
\begin{aligned}
H(\mathbf{k})&=2\lambda_x\sin k_x\sigma_xs_z\tau_z+2\lambda_y\sin k_y\sigma_y\tau_z+(\xi_{\bf k}\sigma_z-\mu)\tau_z+\Delta_0\tau_x+{\bf h\cdot s}\\
\xi_{\bf k}&=\epsilon_0-2t_x\cos k_x-2t_y\cos k_y
\end{aligned}
\end{equation}\label{h1}$$

这里的基矢选择为$\Psi^\dagger=(c^\dagger_{a\uparrow\mathbf{k}},c^\dagger_{b\uparrow\mathbf{k}},c^\dagger_{a\downarrow\mathbf{k}},c^\dagger_{b\downarrow\mathbf{k}},c_{a\downarrow\mathbf{-k}},c_{b\downarrow\mathbf{-k}},-c_{a\uparrow\mathbf{-k}},-c_{b\uparrow\mathbf{-k}})=(C_\mathbf{k}^\dagger,-is_y\sigma_0C_\mathbf{-k})$

在[超导态基矢选择对构建BdG哈密顿量的影响](https://yxli8023.github.io/2021/01/20/BdG-formation2.html)这篇博客中我也提到过,这种基矢的选取虽然从形式上分析是比较方便的,但是在写程序的时候,不仅仅需要考虑哈密顿量中的符号,基矢中的符号也同时需要考虑,根据我自己的习惯,这时候写程序就会变的有些麻烦,还是想选用$\Psi^\dagger=(c^\dagger_{a\uparrow\mathbf{k}},c^\dagger_{b\uparrow\mathbf{k}},c^\dagger_{a\downarrow\mathbf{k}},c^\dagger_{b\downarrow\mathbf{k}},c_{a\downarrow\mathbf{-k}},c_{b\downarrow\mathbf{-k}},c_{a\uparrow\mathbf{-k}},c_{b\uparrow\mathbf{-k}})$这种基矢,这样的话就只需要考虑哈密顿量的具体形式即可,不需要额外再去考虑基矢中的符号.

同样利用[超导态基矢选择对构建BdG哈密顿量的影响](https://yxli8023.github.io/2021/01/20/BdG-formation2.html)这篇文章中的方法,来把BdG形式的构建过程算一下,具体我就不演示这个哈密顿量是怎么做的了,可以参考前面的文章,做法都是完全相同的,因为正常态的基矢选择都是相同的,唯一需要考虑的就是超导的时候,空穴部分的基矢是如何选取的,完成之后,哈密顿量可以改写为

$$\begin{equation}
\begin{aligned}
H^{\mathrm{BdG}}(\mathbf{k})&=(\xi_{\bf k}\sigma_z-\mu +{\bf h\cdot s})\tau_z+2\lambda_x\sin k_x\sigma_xs_z+2\lambda_y\sin k_y\sigma_y\tau_z+\Delta_0s_y\tau_y\\
\xi_{\bf k}&=\epsilon_0-2t_x\cos k_x-2t_y\cos k_y
\end{aligned}
\end{equation}\label{h2}$$

简记$\Gamma_1=\sigma_zs_0\tau_z,\Gamma_2=\sigma_0s_x\tau_z,\Gamma_3=\sigma_xs_z\tau_0,\Gamma_4=\sigma_ys_0\tau_z,\Gamma_5=\sigma_0s_y\tau_y$

# 写程序计算
## Cylinder Geomotery
```fortran
!   https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.227001
!   将基矢写作通常的情况
    module pub
    implicit none
    integer xn,N,hn,kn
    parameter(xn = 50,kn = 50,hn = 8,N = xn*hn)
    real,parameter::pi = 3.1415926535
    complex,parameter::im = (0.0,1.0)  !虚数单位
    complex::Ham(N,N) = 0
!=================================
    real m0   !Driac mass
    real tx,ty
    real lamx,lamy
    real del0
    real h  ! Zeeman field
    complex g1(hn,hn),g2(hn,hn),g3(hn,hn),g4(hn,hn),g5(hn,hn),g6(hn,hn)
!================MKL===============
    integer::lda = N
    integer,parameter::lwmax=2*N+N**2
    real,allocatable::w(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   ! at least 2*N+N**2
    integer lrwork    ! at least 1 + 5*N +2*N**2
    integer liwork   ! at least 3 +5*N
    integer info
    end module pub
!================= PROGRAM START ============================
    program ex01
    use pub
    integer m,l,i
    real ky
    !=====================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1+5*N+2*N**2))
    allocate(iwork(3+5*N))
    !--------------------------
    m0 = 1.0
    tx = 1.0
    ty = 1.0
    lamx = 1.0
    lamy = 1.0
    del0 = 0.4
    h = 0.4
    call band()
    stop 
    end program
!========================================================================
    subroutine band()
    use pub
    real ky
    integer ik
    open(20,file="x-open.dat")
    open(21,file="y-open.dat")
    do ik = -kn,kn
        ky = ik*pi/kn
        call openx(ky)
        write(20,999)ky/pi,(w(i),i = 1,N)
        call openy(ky)
        write(21,999)ky/pi,(w(i),i = 1,N)
    end do
    close(20)
    close(21)
999 format(4001f11.6)
    end subroutine band
!======================== Pauli Matrix driect product============================
    subroutine Pauli()
    use pub
    !----mass term
    g1(1,1) = 1
    g1(2,2) = -1
    g1(3,3) = 1
    g1(4,4) = -1
    g1(5,5) = -1
    g1(6,6) = 1
    g1(7,7) = -1
    g1(8,8) = 1
    !----Zeeman field
    g2(1,3) = 1
    g2(2,4) = 1
    g2(3,1) = 1
    g2(4,2) = 1
    g2(5,7) = -1
    g2(6,8) = -1
    g2(7,5) = -1
    g2(8,6) = -1
    !--- sin(kx)
    g3(1,2) = 1
    g3(2,1) = 1
    g3(3,4) = -1
    g3(4,3) = -1
    g3(5,6) = 1
    g3(6,5) = 1
    g3(7,8) = -1
    g3(8,7) = -1
    !----sin(ky)
    g4(1,2) = -im
    g4(2,1) = im
    g4(3,4) = -im
    g4(4,3) = im
    g4(5,6) = im
    g4(6,5) = -im
    g4(7,8) = im
    g4(8,7) = -im
    !----Delta(k)
    g5(1,7) = -1
    g5(2,8) = -1
    g5(3,5) = 1
    g5(4,6) = 1
    g5(5,3) = 1
    g5(6,4) = 1
    g5(7,1) = -1
    g5(8,2) = -1
    !---mu
    g6(1,1) = 1
    g6(2,2) = 1
    g6(3,3) = 1
    g6(4,4) = 1
    g6(5,5) = -1
    g6(6,6) = -1
    g6(7,7) = -1
    g6(8,8) = -1
    return
    end subroutine Pauli
!==========================================================
    subroutine openx(ky)
    use pub
    real ky
    integer m,l,k
    call Pauli()
    Ham = 0
    !=================== Non-diag Term========================
    do k = 0,xn-1
        if (k == 0) then
            do m = 1,hn
                do l = 1,hn
                    Ham(m,l) = 2*lamy*sin(ky)*g4(m,l) + (m0 - 2*ty*cos(ky))*g1(m,l) + del0*g5(m,l) + h*g2(m,l) + mu*g6(m,l)
                    Ham(m,l + 8) = -im*lamx*g3(m,l) - tx*g1(m,l)
                end do
            end do
        elseif ( k == xn-1 ) then
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) = 2*lamy*sin(ky)*g4(m,l) + (m0 - 2*ty*cos(ky))*g1(m,l) + del0*g5(m,l) + h*g2(m,l) + mu*g6(m,l)
                    Ham(k*hn + m,k*hn + l - hn) = im*lamx*g3(m,l) - tx*g1(m,l)
                end do
            end do
        else
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) = 2*lamy*sin(ky)*g4(m,l) + (m0 - 2*ty*cos(ky))*g1(m,l) + del0*g5(m,l) + h*g2(m,l) + mu*g6(m,l)
                    Ham(k*hn + m,k*hn + l + hn) = -im*lamx*g3(m,l) - tx*g1(m,l)
                    Ham(k*hn + m,k*hn + l - hn) = im*lamx*g3(m,l) - tx*g1(m,l)
                end do
            end do
        end if
    end do
    call ishermitian()
    call eigsol()
    end subroutine openx
!==========================================================
    subroutine openy(kx)
    use pub
    real kx
    integer m,l,k
    call Pauli()
    Ham = 0
    !=================== Non-diag Term========================
    do k = 0,xn-1
        if (k == 0) then
            do m = 1,hn
                do l = 1,hn
                    Ham(m,l) = 2*lamx*sin(kx)*g3(m,l) + (m0 - 2*tx*cos(kx))*g1(m,l) + del0*g5(m,l) + h*g2(m,l) + mu*g6(m,l)
                    Ham(m,l + 8) = -im*lamy*g4(m,l) - ty*g1(m,l)
                end do
            end do
        elseif ( k == xn-1 ) then
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) = 2*lamx*sin(kx)*g3(m,l) + (m0 - 2*tx*cos(kx))*g1(m,l) + del0*g5(m,l) + h*g2(m,l) + mu*g6(m,l)
                    Ham(k*hn + m,k*hn + l - hn) = im*lamy*g4(m,l) - ty*g1(m,l)
                end do
            end do
        else
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) = 2*lamx*sin(kx)*g3(m,l) + (m0 - 2*tx*cos(kx))*g1(m,l) + del0*g5(m,l) + h*g2(m,l) + mu*g6(m,l)
                    Ham(k*hn + m,k*hn + l + hn) = -im*lamy*g4(m,l) - tx*g1(m,l)
                    Ham(k*hn + m,k*hn + l - hn) = im*lamy*g4(m,l) - tx*g1(m,l)
                end do
            end do
        end if
    end do
    call ishermitian()
    call eigsol()
    end subroutine openy
!============================================================
    subroutine ishermitian()
    use pub
    integer i,j
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                open(160,file = 'hermitian.dat')
                ccc = ccc +1
                write(160,*)i,j
                write(160,*)Ham(i,j)
                write(160,*)Ham(j,i)
                write(160,*)"===================="
                stop
            end if
        end do
    end do
    close(160)
    return
    end subroutine ishermitian
!================= 矩阵本征值求解 ==============
    subroutine eigSol()
    use pub
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(110,file="mes.txt",status="unknown")
        write(110,*)'The algorithm failed to compute eigenvalues.'
        close(110)
    end if
    open(120,file="eigval.dat",status="unknown")
    do m = 1,N
        write(120,*)m,w(m)
    end do
    close(120)
    return
    end subroutine eigSol
```
这个程序是用来计算边界态的,即一个方向是开边界,另外一个方向是周期边界条件,我已验证过,相同的参数下和文章中结果是相同.
## Real-space
```fortran
!   https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.227001
!   将基矢写作通常的情况
    module pub
    implicit none
    integer xn,N,yn,hn,kn,len2
    parameter(xn = 36,yn = xn,kn = 50,hn = 8,len2 = xn*yn,N = len2*hn)
    real eps,pi
    parameter (pi = 3.1415926535, eps = 1e-5)
    complex,parameter::im = (0.0,1.0)  !虚数单位
    complex::Ham(N,N) = 0
    !-----------------------------------------
    real m0,tx,ty,lamx,lamy,del0,h
    integer bry(4,len2)  ! boundary index
    !------------------------------------------
    integer::lda = N
    integer,parameter::lwmax=2*N+N**2
    real,allocatable::w(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   ! at least 2*N+N**2
    integer lrwork    ! at least 1 + 5*N +2*N**2
    integer liwork   ! at least 3 +5*N
    integer info
    end module pub
!================= PROGRAM START ============================
    program ex01
    use pub
    integer m,l,i
    !=====================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1+5*N+2*N**2))
    allocate(iwork(3+5*N))
    !--------------------------
    m0 = 1.0
    tx = 1.0
    ty = 1.0
    lamx = 1.0
    lamy = 1.0
    del0 = 0.4
    h = 0.8
    call matset2()
    call ldos()
    stop 
    end program
!===========================Local Density of State=============================
    subroutine ldos()
    use pub
    integer ix,iy
    real,external::wave
    open(12,file="wavenorm.dat") 
    do iy = 1,yn
        do ix = 1,xn
            write(12,*)ix,iy,wave(ix,iy)
        end do
    end do
    close(12)
    return
    end subroutine ldos
!=======================================================================
    real function wave(x,y)
    use pub
    integer x,y,ind
    ind = (y-1)*xn + x
    wave = 0
    do i2 = 0,7
        do zo = N/2 - 3,N/2 + 4
            wave = wave + abs(Ham(ind + len2*i2,zo))**2
        end do
    end do
    end function wave
!==========================================================
    subroutine boundary()
    use pub
    integer ix,iy,i
    bry = 0
    do iy = 1,yn
        do ix = 1,xn
            i = (iy - 1)*xn + ix
            bry(1,i) = i + 1
            if(ix == xn)bry(1,i) = bry(1,i) - xn
            bry(2,i) = i - 1
            if(ix == 1)bry(2,i) = bry(2,i) + xn
            bry(3,i) = i + yn
            if(iy == yn)bry(3,i) = bry(3,i) - len2
            bry(4,i) = i - yn
            if(iy == 1)bry(4,i) = bry(4,i) + len2
        end do
    end do
    return
    end subroutine boundary
!=====================================================================
    subroutine matset2()
    use pub
    integer ix,iy,i
    call boundary()
    Ham = 0
    do iy = 1,yn
        do ix = 1,xn
            i = (iy - 1)*xn + ix
            !-------------------------------
            !  Mass term
            !(1,1)
            ham(i,i) = m0 + mu
            if(ix.ne.xn)ham(i,bry(1,i)) = tx
            if(ix.ne.1)ham(i,bry(2,i)) = tx
            if(iy.ne.yn)ham(i,bry(3,i)) = ty
            if(iy.ne.1)ham(i,bry(4,i)) = ty
            !(2,2)
            ham(len2 + i,len2 + i) = -m0 + mu
            if(ix.ne.xn)ham(len2 + i,len2 + bry(1,i)) = -tx
            if(ix.ne.1)ham(len2 + i,len2 + bry(2,i)) = -tx
            if(iy.ne.yn)ham(len2 + i,len2 + bry(3,i)) = -ty
            if(iy.ne.1)ham(len2 + i,len2 + bry(4,i)) = -ty
            !(3,3)
            ham(len2*2 + i,len2*2 + i) = m0 + mu
            if(ix.ne.xn)ham(len2*2 + i,len2*2 + bry(1,i)) = tx
            if(ix.ne.1)ham(len2*2 + i,len2*2 + bry(2,i)) = tx
            if(iy.ne.yn)ham(len2*2 + i,len2*2 + bry(3,i)) = ty
            if(iy.ne.1)ham(len2*2 + i,len2*2 + bry(4,i)) = ty
            !(4,4)
            ham(len2*3 + i,len2*3 + i) = -m0 + mu
            if(ix.ne.xn)ham(len2*3 + i,len2*3 + bry(1,i)) = -tx
            if(ix.ne.1)ham(len2*3 + i,len2*3 + bry(2,i)) = -tx
            if(iy.ne.yn)ham(len2*3 + i,len2*3 + bry(3,i)) = -ty
            if(iy.ne.1)ham(len2*3 + i,len2*3 + bry(4,i)) = -ty
            !(5,5)
            ham(len2*4 + i,len2*4 + i) = -m0 - mu
            if(ix.ne.xn)ham(len2*4 + i,len2*4 + bry(1,i)) = -tx
            if(ix.ne.1)ham(len2*4 + i,len2*4 + bry(2,i)) = -tx
            if(iy.ne.yn)ham(len2*4 + i,len2*4 + bry(3,i)) = -ty
            if(iy.ne.1)ham(len2*4 + i,len2*4 + bry(4,i)) = -ty
            !(6,6)
            ham(len2*5 + i,len2*5 + i) = m0 - mu
            if(ix.ne.xn)ham(len2*5 + i,len2*5 + bry(1,i)) = tx
            if(ix.ne.1)ham(len2*5 + i,len2*5 + bry(2,i)) = tx
            if(iy.ne.yn)ham(len2*5 + i,len2*5 + bry(3,i)) = ty
            if(iy.ne.1)ham(len2*5 + i,len2*5 + bry(4,i)) = ty
            !(7,7)
            ham(len2*6 + i,len2*6 + i) = -m0 - mu
            if(ix.ne.xn)ham(len2*6 + i,len2*6 + bry(1,i)) = -tx
            if(ix.ne.1)ham(len2*6 + i,len2*6 + bry(2,i)) = -tx
            if(iy.ne.yn)ham(len2*6 + i,len2*6 + bry(3,i)) = -ty
            if(iy.ne.1)ham(len2*6 + i,len2*6 + bry(4,i)) = -ty
            !(8,8)
            ham(len2*7 + i,len2*7 + i) = m0 - mu
            if(ix.ne.xn)ham(len2*7 + i,len2*7 + bry(1,i)) = tx
            if(ix.ne.1)ham(len2*7 + i,len2*7 + bry(2,i)) = tx
            if(iy.ne.yn)ham(len2*7 + i,len2*7 + bry(3,i)) = ty
            if(iy.ne.1)ham(len2*7 + i,len2*7 + bry(4,i)) = ty
            !----------------------------------------------
            !   SOC
            !(1,2)
            if(ix.ne.xn)ham(i,len2 + bry(1,i)) = lamx/im
            if(ix.ne.1)ham(i,len2 + bry(2,i)) = -lamx/im
            if(iy.ne.yn)ham(i,len2 + bry(3,i)) = -im*lamx/im
            if(iy.ne.1)ham(i,len2 + bry(4,i)) = im*lamx/im
            !(2,1)
            if(ix.ne.xn)ham(len2 + i,bry(1,i)) = lamx/im
            if(ix.ne.1)ham(len2 + i,bry(2,i)) = -lamx/im
            if(iy.ne.yn)ham(len2 + i,bry(3,i)) = im*lamx/im
            if(iy.ne.1)ham(len2 + i,bry(4,i)) = -im*lamx/im
            !(3,4)
            if(ix.ne.xn)ham(len2*2 + i,len2*3 + bry(1,i)) = -lamx/im
            if(ix.ne.1)ham(len2*2 + i,len2*3 + bry(2,i)) = lamx/im
            if(iy.ne.yn)ham(len2*2 + i,len2*3 + bry(3,i)) = -im*lamx/im
            if(iy.ne.1)ham(len2*2 + i,len2*3 + bry(4,i)) = im*lamx/im
            !(4,3)
            if(ix.ne.xn)ham(len2*3 + i,len2*2 + bry(1,i)) = -lamx/im
            if(ix.ne.1)ham(len2*3 + i,len2*2 + bry(2,i)) = lamx/im
            if(iy.ne.yn)ham(len2*3 + i,len2*2 + bry(3,i)) = im*lamx/im
            if(iy.ne.1)ham(len2*3 + i,len2*2 + bry(4,i)) = -im*lamx/im
            !(5,6)
            if(ix.ne.xn)ham(len2*4 + i,len2*5 + bry(1,i)) = lamx/im
            if(ix.ne.1)ham(len2*4 + i,len2*5 + bry(2,i)) = -lamx/im
            if(iy.ne.yn)ham(len2*4 + i,len2*5 + bry(3,i)) = im*lamx/im
            if(iy.ne.1)ham(len2*4 + i,len2*5 + bry(4,i)) = -im*lamx/im
            !(6,5)
            if(ix.ne.xn)ham(len2*5 + i,len2*4 + bry(1,i)) = lamx/im
            if(ix.ne.1)ham(len2*5 + i,len2*4 + bry(2,i)) = -lamx/im
            if(iy.ne.yn)ham(len2*5 + i,len2*4 + bry(3,i)) = -im*lamx/im
            if(iy.ne.1)ham(len2*5 + i,len2*4 + bry(4,i)) = im*lamx/im
            !(7,8)
            if(ix.ne.xn)ham(len2*6 + i,len2*7 + bry(1,i)) = -lamx/im
            if(ix.ne.1)ham(len2*6 + i,len2*7 + bry(2,i)) = lamx/im
            if(iy.ne.yn)ham(len2*6 + i,len2*7 + bry(3,i)) = im*lamx/im
            if(iy.ne.1)ham(len2*6 + i,len2*7 + bry(4,i)) = -im*lamx/im
            !(8,7)
            if(ix.ne.xn)ham(len2*7 + i,len2*6 + bry(1,i)) = -lamx/im
            if(ix.ne.1)ham(len2*7 + i,len2*6 + bry(2,i)) = lamx/im
            if(iy.ne.yn)ham(len2*7 + i,len2*6 + bry(3,i)) = -im*lamx/im
            if(iy.ne.1)ham(len2*7 + i,len2*6 + bry(4,i)) = im*lamx/im
            !---------------------------------------------------
            !  SC
            !(1,7)
            ham(i,len2*6 + i) = -del0
            !(2,8)
            ham(len2 + i,len2*7 + i) = -del0
            !(3,5)
            ham(len2*2 + i,len2*4 + i) = del0
            !(4,6)
            ham(len2*3 + i,len2*5 + i) = del0
            !(5,3)
            ham(len2*4 + i,len2*2 + i) = del0
            !(6,4)
            ham(len2*5 + i,len2*3 + i) = del0
            !(7,1)
            ham(len2*6 + i,i) = -del0
            !(8,2)
            ham(len2*7 + i,len2 + i) = -del0
            !---------------------------------------------------
            ! Zeeman field
            !(1,3)
            ham(i,len2*2 + i) = h
            !(2,4)
            ham(len2 + i,len2*3 + i) = h
            !(3,1)
            ham(len2*2 + i,i) = h
            !(4,2)
            ham(len2*3 + i,len2 + i) = h
            !(5,7)
            ham(len2*4 + i,len2*6 + i) = -h
            !(6,8)
            ham(len2*5 + i,len2*7 + i) = -h
            !(7,5)
            ham(len2*6 + i,len2*4 + i) = -h
            !(8,6)
            ham(len2*7 + i,len2*5 + i) = -h
        end do
    end do
    !-------------------------
    call ishermitian()
    call eigsol()
    end subroutine matset2
!============================================================
    subroutine ishermitian()
    use pub
    integer i,j
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                open(160,file = 'hermitian.dat')
                ccc = ccc +1
                write(160,*)i,j
                write(160,*)Ham(i,j)
                write(160,*)Ham(j,i)
                write(160,*)"===================="
                write(*,*)"Ham isn't Hermitian"
                stop
            end if
        end do
    end do
    close(160)
    return
    end subroutine ishermitian
!================= 矩阵本征值求解 ==============
    subroutine eigSol()
    use pub
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(110,file="mes.dat",status="unknown")
        write(110,*)'The algorithm failed to compute eigenvalues.'
        close(110)
    end if
    open(120,file="eigval.dat",status="unknown")
    do m = 1,N
        write(120,*)m,w(m)
    end do
    close(120)
    return
    end subroutine eigSol
```
这个程序是用来计算实空间波函数分布的,即两个方向都是开边界条件,我已验证过,相同的参数下和文章中结果是相同.

所有程序的编译命令为
> ifort -mkl file.f90 -o a.out
> ./a.out &   ! 后台执行程序


# 参考
1.[Bogoliubov变换与Majorana表示](https://zhuanlan.zhihu.com/p/59445571)

2.[Bogoliubov-de Gennes Method and Its Applications](https://link.springer.com/book/10.1007/978-3-319-31314-6)

3.[High-Temperature Majorana Corner States](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.186801)

4.[Majorana Corner Modes in a High-Temperature Platform](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.096803)

5.[n-Plane Zeeman-Field-Induced Majorana Corner and Hinge Modes in an s-Wave Superconductor Heterostructure](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.227001)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)