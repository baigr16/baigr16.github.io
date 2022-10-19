---
title: Fortran计算选择区间内本征值
tags: Fortran 
layout: article
license: true
toc: true
key: a20220826
pageview: true
# cover: /assets/images/GroupTheory/cube_symmetry.jpg
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
在进行数值计算的时候，对于较大的矩阵，有时候并不完全需要其所有的本征值，而只需要一定范围内的本征值即可，这里就整理一下如何利用`Fortran`来得到选定区间内的本征值
{:.info}
<!--more-->
这里其实就是使用`CHEEVX`这个函数进行厄米矩阵对角化过程，具体使用实例可以参考[这里](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top/least-squares-and-eigenvalue-problems/symmetric-eigenproblems/heevx-function/cheevx-example/cheevx-example-fortran.html)，该函数的帮助文档在[这里](https://netlib.org/lapack/explore-html/d9/de3/group__complex_h_eeigen_ga9f7c713a0119e777afe726e54feb6ef7.html)。
# 代码
```fortran
!   Author:YuXuan-Li
!   Email:yxli406@gmail.com
!   Ref:https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.107.097001
    module pub
    implicit none
    integer xn,yn,kn,N,len2,en
    parameter(xn = 30, yn = xn, kn = 100,N = xn*yn*8,len2 = xn*yn,en = 200)
    real,parameter::pi = 3.1415926535
    complex,parameter::im = (0.,1.0)  
    integer bry(4,len2)
    complex Ham(N,N) 
    real tx,ty,tz,m0,mx,my,mz,del0,mu,chi,delx,dely
    !-------------------------------------
    integer::lda = N
    integer::ldz = N
    integer,parameter::lwmax = 2*N + N**2
    integer INFO, LWORK, IL,IU, M
    real ABSTOL,VL,VU
    real,allocatable::w(:)
    integer IFAIL(N)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    complex,allocatable:: Z(:,:)
    end module pub
!============================================================
    program sol
    use pub
    real t1,t2
    !======================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1+5*N+2*N**2))
    allocate(iwork(3+5*N))
    allocate(Z(LDZ,N))
    !----------------------
    tx = 0.5
    ty = 0.5
    tz = 0.5
    m0 = 2.5
    mx = -1.0
    my = -1.0
    mz = -1.0
    delx = 0
    dely = 0
    del0 = 0.1
    mu = 0.0
    chi = 1.0
    ! call matrixSet(0.0)
    ! call eigSol()
    !-----------------------
    ABSTOL = -1.0 !Negative ABSTOL means using the default value
    !Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
    VL = -0.1
    VU = 0.1
    !-----------------------
    call cpu_time(t1)
    call matrixSet(0.0)
    call eigSol2()
    call cpu_time(t2)
    write(*,*)t2 - t1
    call cpu_time(t1)
    call matrixSet(0.0)
    call eigSol1()
    call cpu_time(t2)
    write(*,*)t2 - t1
    ! call muchange()
    stop
    end program sol
!============================================================
    subroutine muchange()
    use pub
    integer m1
    open(33,file="mu.dat")
    do m1 = 1,en
        mu = m1*1.5/en
        call matrixSet(0.0)
        call eigSol2()
        write(33,999)mu,(w(i),i=3550,3650)
    end do
    close(33)
999 format(110f11.6)
    return
    end subroutine muchange
!============================================================
    subroutine  matrixSet(kz)
    !   Bais: a-up b-up a-down b-down
    use pub
    integer ix,iy,i
    real dem,num,kz
    complex theta
    call boundary()
    Ham = 0
    do iy = 1,yn
        do ix = 1,xn
            i = (iy-1)*xn + ix
            !============= Kinetic
            Ham(i,i) = m0 - mu + mz*cos(kz)
            Ham(len2 + i,len2 + i) = -m0 - mu - mz*cos(kz)
            Ham(len2*2 + i,len2*2 + i) = m0 - mu + mz*cos(kz)
            Ham(len2*3 + i,len2*3 + i) = -m0 - mu - mz*cos(kz)
            Ham(len2*4 + i,len2*4 + i) = -m0 + mu - mz*cos(kz)
            Ham(len2*5 + i,len2*5 + i) = m0 + mu + mz*cos(kz)
            Ham(len2*6 + i,len2*6 + i) = -m0 + mu - mz*cos(kz)
            Ham(len2*7 + i,len2*7 + i) = m0 + mu + mz*cos(kz)
            !----
            !(1,1)
            if(ix.ne.xn)Ham(i,bry(1,i)) = mx/2.0
            if(ix.ne.1)Ham(i,bry(2,i)) = mx/2.0
            if(iy.ne.yn)Ham(i,bry(3,i)) = my/2.0
            if(iy.ne.1)Ham(i,bry(4,i)) = my/2.0
            !(2,2)
            if(ix.ne.xn)Ham(len2 + i,len2 + bry(1,i)) = -mx/2.0
            if(ix.ne.1)Ham(len2 + i,len2 + bry(2,i)) = -mx/2.0
            if(iy.ne.yn)Ham(len2 + i,len2 + bry(3,i)) = -my/2.0
            if(iy.ne.1)Ham(len2 + i,len2 + bry(4,i)) = -my/2.0
            !(3,3)
            if(ix.ne.xn)Ham(len2*2 + i,len2*2 + bry(1,i)) = mx/2.0
            if(ix.ne.1)Ham(len2*2 + i,len2*2 + bry(2,i)) = mx/2.0
            if(iy.ne.yn)Ham(len2*2 + i,len2*2 + bry(3,i)) = my/2.0
            if(iy.ne.1)Ham(len2*2 + i,len2*2 + bry(4,i)) = my/2.0
            !(4,4)
            if(ix.ne.xn)Ham(len2 *3 + i,len2*3 + bry(1,i)) = -mx/2.0
            if(ix.ne.1)Ham(len2 *3 + i,len2*3 + bry(2,i)) = -mx/2.0
            if(iy.ne.yn)Ham(len2 *3 + i,len2*3 + bry(3,i)) = -my/2.0
            if(iy.ne.1)Ham(len2 *3 + i,len2*3 + bry(4,i)) = -my/2.0
            !(5,5)
            if(ix.ne.xn)Ham(len2*4 + i,len2*4 + bry(1,i)) = -mx/2.0
            if(ix.ne.1)Ham(len2*4 + i,len2*4 + bry(2,i)) = -mx/2.0
            if(iy.ne.yn)Ham(len2*4 + i,len2*4 + bry(3,i)) = -my/2.0
            if(iy.ne.1)Ham(len2*4 + i,len2*4 + bry(4,i)) = -my/2.0
            !(6,6)
            if(ix.ne.xn)Ham(len2*5 + i,len2*5 + bry(1,i)) = mx/2.0
            if(ix.ne.1)Ham(len2*5 + i,len2*5 + bry(2,i)) = mx/2.0
            if(iy.ne.yn)Ham(len2*5 + i,len2*5 + bry(3,i)) = my/2.0
            if(iy.ne.1)Ham(len2*5 + i,len2*5 + bry(4,i)) = my/2.0
            !(7,7)
            if(ix.ne.xn)Ham(len2*6 + i,len2*6 + bry(1,i)) = -mx/2.0
            if(ix.ne.1)Ham(len2*6 + i,len2*6 + bry(2,i)) = -mx/2.0
            if(iy.ne.yn)Ham(len2*6 + i,len2*6 + bry(3,i)) = -my/2.0
            if(iy.ne.1)Ham(len2*6 + i,len2*6 + bry(4,i)) = -my/2.0
            !(8,8)
            if(ix.ne.xn)Ham(len2*7 + i,len2*7 + bry(1,i)) = mx/2.0
            if(ix.ne.1)Ham(len2*7 + i,len2*7 + bry(2,i)) = mx/2.0
            if(iy.ne.yn)Ham(len2*7 + i,len2*7 + bry(3,i)) = my/2.0
            if(iy.ne.1)Ham(len2*7 + i,len2*7 + bry(4,i)) = my/2.0
            !===========SOC
            !(1,4)
            if(ix.ne.xn)Ham(i,len2*3 + bry(1,i)) = -im*tx
            if(ix.ne.1)Ham(i,len2*3 + bry(2,i)) = im*tx
            if(iy.ne.yn)Ham(i,len2*3 + bry(3,i)) = -ty
            if(iy.ne.1)Ham(i,len2*3 + bry(4,i)) = ty
            !(2,3)
            if(ix.ne.xn)Ham(len2 + i,len2*2 + bry(1,i)) = -im*tx
            if(ix.ne.1)Ham(len2 + i,len2*2 + bry(2,i)) = im*tx
            if(iy.ne.yn)Ham(len2 + i,len2*2 + bry(3,i)) = ty
            if(iy.ne.1)Ham(len2 + i,len2*2 + bry(4,i)) = -ty
            !(3,2)
            if(ix.ne.xn)Ham(len2*2 + i,len2 + bry(1,i)) = -im*tx
            if(ix.ne.1)Ham(len2*2 + i,len2 + bry(2,i)) = im*tx
            if(iy.ne.yn)Ham(len2*2 + i,len2 + bry(3,i)) = -ty
            if(iy.ne.1)Ham(len2*2 + i,len2 + bry(4,i)) = ty
            !(4,1)
            if(ix.ne.xn)Ham(len2*3 + i,bry(1,i)) = -im*tx
            if(ix.ne.1)Ham(len2*3 + i,bry(2,i)) = im*tx
            if(iy.ne.yn)Ham(len2*3 + i,bry(3,i)) = ty
            if(iy.ne.1)Ham(len2*3 + i,bry(4,i)) = -ty
            !(5,8)
            Ham(len2*4 + i,len2*7 + bry(1,i)) = -conjg(Ham(i,len2*3 + bry(1,i)))
            Ham(len2*4 + i,len2*7 + bry(2,i)) = -conjg(Ham(i,len2*3 + bry(2,i)))
            Ham(len2*4 + i,len2*7 + bry(3,i)) = -conjg(Ham(i,len2*3 + bry(3,i)))
            Ham(len2*4 + i,len2*7 + bry(4,i)) = -conjg(Ham(i,len2*3 + bry(4,i)))
            !(6,7)
            Ham(len2*5 + i,len2*6 + bry(1,i)) = -conjg(Ham(len2 + i,len2*2 + bry(1,i)))
            Ham(len2*5 + i,len2*6 + bry(2,i)) = -conjg(Ham(len2 + i,len2*2 + bry(2,i)))
            Ham(len2*5 + i,len2*6 + bry(3,i)) = -conjg(Ham(len2 + i,len2*2 + bry(3,i)))
            Ham(len2*5 + i,len2*6 + bry(4,i)) = -conjg(Ham(len2 + i,len2*2 + bry(4,i)))
            !(7,6)
            Ham(len2*6 + i,len2*5 + bry(1,i)) = -conjg(Ham(len2*2 + i,len2 + bry(1,i)))
            Ham(len2*6 + i,len2*5 + bry(2,i)) = -conjg(Ham(len2*2 + i,len2 + bry(2,i)))
            Ham(len2*6 + i,len2*5 + bry(3,i)) = -conjg(Ham(len2*2 + i,len2 + bry(3,i)))
            Ham(len2*6 + i,len2*5 + bry(4,i)) = -conjg(Ham(len2*2 + i,len2 + bry(4,i)))
            !(8,5)
            Ham(len2*7 + i,len2*4 + bry(1,i)) = -conjg(Ham(len2*3 + i,bry(1,i)))
            Ham(len2*7 + i,len2*4 + bry(2,i)) = -conjg(Ham(len2*3 + i,bry(2,i)))
            Ham(len2*7 + i,len2*4 + bry(3,i)) = -conjg(Ham(len2*3 + i,bry(3,i)))
            Ham(len2*7 + i,len2*4 + bry(4,i)) = -conjg(Ham(len2*3 + i,bry(4,i)))
            !===============================================================
            Ham(i,len2 + i) = tz*sin(kz)
            Ham(len2 + i,i)  = tz*sin(kz)
            Ham(len2*2 + i,len2*3 + i) = -tz*sin(kz)
            Ham(len2*3 + i,len2*2 + i) = -tz*sin(kz)

            Ham(len2*4 + i,len2*5 + i) = -tz*sin(kz)
            Ham(len2*5 + i,len2*4 + i) = -tz*sin(kz)
            Ham(len2*6 + i,len2*7 + i) = tz*sin(kz)
            Ham(len2*7 + i,len2*6 + i) = tz*sin(kz)

            !========SC Pairing
            dem = (ix - 15)*1.0
            num = (iy - 15)*1.0
            theta = tanh(sqrt(dem**2 + num**2)/chi)*exp(im*atan(num/(dem + 0.00001)))
            ! write(33,*)ix,iy,real(theta),aimag(theta)
            Ham(i,len2*6 + i) = del0*exp(im*theta)
            Ham(len2 + i,len2*7 + i) = del0*exp(im*theta)
            Ham(len2*2 + i,len2*4 + i) = -del0*exp(im*theta)
            Ham(len2*3 + i,len2*5 + i) = -del0*exp(im*theta)

            Ham(len2*4 + i,len2*2 + i) = conjg(Ham(len2*2 + i,len2*4 + i))
            Ham(len2*5 + i,len2*3 + i) = conjg(Ham(len2*3 + i,len2*5 + i))
            Ham(len2*6 + i,i) = conjg(Ham(i,len2*6 + i))
            Ham(len2*7 + i,len2 + i) = conjg(Ham(len2 + i,len2*7 + i))
        end do
    end do
    !------------------------------
    call isHermitian()
    return
    end subroutine matrixSet
!============================================================
    subroutine boundary()
    use pub
    integer ix,iy,i
    do iy = 1,yn
        do ix = 1,xn
            i = (iy-1)*xn + ix
            bry(1,i) = i + 1
            if(ix==xn)bry(1,i) = bry(1,i) - xn
            bry(2,i) = i - 1
            if(ix==1)bry(2,i) = bry(2,i) + xn
            bry(3,i) = i + xn
            if(iy==yn)bry(3,i) = bry(3,i) - len2
            bry(4,i) = i - xn
            if(iy==1)bry(4,i) = bry(4,i) + len2
        end do
    end do
    return
    end subroutine boundary
!============================================================
    subroutine isHermitian()
    use pub
    integer i,j
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                open(16,file = 'hermitian.dat')
                write(16,*)i,j
                write(16,*)Ham(i,j)
                write(16,*)Ham(j,i)
                write(16,*)"===================="
                write(*,*)"Hamiltonian is not Hermitian"
                stop
            end if
        end do
    end do
    close(16)
    return
    end subroutine isHermitian
!================= 矩阵本征值求解 ==============
    subroutine eigSol1()
    use pub
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork &
        ,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','Upper',N,Ham,lda,w,work,lwork &
        ,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(11,file="mes.dat",status="unknown")
        write(11,*)'The algorithm failed to compute eigenvalues.'
        close(11)
    end if
    open(120,file="eigval1.dat")
    do i0 = 1,N
        write(120,*)i0,w(i0)
    end do
    close(120)
    return
    end subroutine eigSol1
!==================================================
    subroutine eigSol2()
    use pub
    integer i0
    lwork = -1
    liwork = -1
    lrwork = -1
    LWORK = -1
    CALL CHEEVX( 'Vectors', 'Values', 'Lower', N, Ham, LDA, VL, VU, IL, &
                IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK,IFAIL, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
!     Solve eigenproblem.
    CALL CHEEVX( 'Vectors', 'Values', 'Lower', N, Ham, LDA, VL, VU, IL,&
                IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO )
    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        STOP
    END IF
    open(120,file="eigval2.dat")
    do i0 = 1,N
        write(120,*)i0,w(i0)
    end do
    close(120)
    return
    end subroutine eigSol2
```
通过求解这个$30\times 30\times 8$的矩阵，发现如果只是求解少数的本征值，计算速度会得到显著的提升。
# 参考
- [CHEEVD函数文档](https://netlib.org/lapack/explore-html/d9/de3/group__complex_h_eeigen_ga9f7c713a0119e777afe726e54feb6ef7.html)
- [CHEEVD实例](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top/least-squares-and-eigenvalue-problems/symmetric-eigenproblems/heevx-function/cheevx-example/cheevx-example-fortran.html)
