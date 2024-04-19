---
title: Fortran结合MPI并行
tags:  Fortran Code MPI 
layout: article
license: true
toc: true
key: a20240419c
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
该Bolg整理在`Fortran`用`MPI`编写并行程序，用在科学计算中速度提升肉眼可见，
{:.info}
<!--more-->
# 代码解析
先给一个完整的代码，是用来计算自旋极化率的，实际上在[Fortran结合MPI并行计算自旋极化率](https://yxli8023.github.io/2024/04/19/Frotran-MPI.html)中已经出现过了，这里就是解读一下其中的并行部分是如何用`MPI`写的。

```fortran
module param
    implicit none
    integer kn,hn
    parameter(kn = 64,hn = 4)
    real,parameter::pi = 3.1415926535897
    real,parameter::omega = 0.00
    complex,parameter::im = (0.,1.) !Imagine unit
    complex Ham(hn,hn),Umat(hn,hn),ones(hn,hn) ! Hamiltonian and interaction Matrix
    real t0,t1x,t1z,t2x,t2z,t3xz,t4xz,tvx,tvz,ex,ez  ! 哈密顿量参数
    real U0,J0  ! 相互作用参数
    real valmesh(2,-kn:kn-1,-kn:kn-1,hn)
    complex vecmesh(2,-kn:kn-1,-kn:kn-1,hn,hn)
    complex chi0(-kn:kn-1,-kn:kn-1,hn,hn),chi(-kn:kn-1,-kn:kn-1,hn,hn)
    ! LAPACK PACKAGE PARAM
    integer::lda = hn
    integer,parameter::lwmax = 2*hn+hn**2
    real,allocatable::val(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   ! at least 2*hn+N**2
    integer lrwork    ! at least 1 + 5*hn +2*hn**2
    integer liwork   ! at least 3 +5*hn
    integer info
end module param
!================================================================================================================================================================================================
program main
    use param
    use mpi 
    !------------------------------------
    complex temp1(hn,hn),temp2(hn,hn),temp3
    integer iky,ikx
    real qx,qy,local_rechi(-kn:kn-1,-kn:kn-1,4),rechi(-kn:kn-1,-kn:kn-1,4)
    !---------------------------------
    integer  numcore,indcore,ierr
    character(len = 20)::filename,char1,char2
    !---------------------------------
    external::cheevd
    allocate(val(hn))
    allocate(work(lwmax))
    allocate(rwork(1+5*hn+2*hn**2))
    allocate(iwork(3+5*hn))
    !---------------------------------
    call MPI_INIT(ierr)     ! 初始化进程
    call MPI_COMM_RANK(MPI_COMM_WORLD, indcore, ierr) ! 得到本进程在通信空间中的rank值,即在组中的逻辑编号(该 rank值为0到p-1间的整数,相当于进程的ID。)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numcore, ierr) !获得进程个数 size, 这里用变量p保存
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    nki = floor(indcore * (2.0 * kn)/numcore) - kn
    nkf = floor((indcore + 1) * (2.0 * kn)/numcore) - kn - 1 
    ! 遍历BZ求极化率
    do iky = nki,nkf
        qy = pi * iky/kn
        do ikx = -kn,kn - 1
            qx = pi * ikx/kn
            call chi0cal(qx,qy,chi0(iky,ikx,:,:))  ! 得到裸极化率
            call inv(ones - matmul(chi0(iky,ikx,:,:),Umat),temp2)  ! 矩阵求逆
            chi(iky,ikx,:,:) = matmul(temp2,chi0(iky,ikx,:,:))
            temp3 = sum(chi(iky,ikx,:,:))
            local_rechi(iky,ikx,1) = qx
            local_rechi(iky,ikx,2) = qy
            local_rechi(iky,ikx,3) = real(temp3)
            local_rechi(iky,ikx,4) = aimag(temp3)
        end do
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! call MPI_Gather(local_rechi, 2 * kn, MPI_COMPLEX, rechi, 2 * kn, MPI_COMPLEX, 0, MPI_COMM_WORLD,ierr)
    call MPI_Reduce(local_rechi, rechi, (2 * kn)**2 * 4, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
    if (indcore .eq. 0) then
        char1 = "fortran-chi-"
        write(char2,"(I3.3)")2 * kn
        filename = trim(char1)//trim(char2)
        char1 = ".dat"
        filename = trim(filename)//trim(char1)
        open(12,file = filename)
        do iky = -kn,kn - 1
            do ikx = -kn,kn - 1
                write(12,"(4F8.3)")rechi(iky,ikx,1),rechi(iky,ikx,2),rechi(iky,ikx,3),rechi(iky,ikx,4)
            end do
        end do
        close(12)
    end if
    call MPI_Finalize(ierr)
    stop
end program main
!================================================================================================================================================================================================
subroutine matset(kx,ky)
    ! 矩阵赋值
    use param
    real kx,ky
    integer k0
    t0 = 1.0
    t1x = -0.483 * t0
    t1z = -0.110 * t0
    t2x = 0.069 * t0
    t2z = -0.017 * t0
    t3xz = 0.239 * t0
    t4xz = -0.034 * t0
    tvx = 0.005 * t0
    tvz = -0.635 * t0
    ex = 0.776 * t0
    ez = 0.409 * t0
    
    Ham = 0.0
    Ham(1, 1) = 2 * t1x * (cos(kx) + cos(ky)) + 4 * t2x * cos(kx) * cos(ky) + ex
    Ham(2, 2) = 2 * t1z * (cos(kx) + cos(ky)) + 4 * t2z * cos(kx) * cos(ky) + ez
    Ham(1, 2) = 2 * t3xz * (cos(kx) - cos(ky))
    Ham(2, 1) = 2 * t3xz * (cos(kx) - cos(ky))
    Ham(3, 3) = 2 * t1x * (cos(kx) + cos(ky)) + 4 * t2x * cos(kx) * cos(ky) + ex
    Ham(4, 4) = 2 * t1z * (cos(kx) + cos(ky)) + 4 * t2z * cos(kx) * cos(ky) + ez
    Ham(3, 4) = 2 * t3xz * (cos(kx) - cos(ky))
    Ham(4, 3) = 2 * t3xz * (cos(kx) - cos(ky))
    Ham(1, 3) = tvx
    Ham(1, 4) = 2 * t4xz * (cos(kx) - cos(ky))
    Ham(2, 3) = 2 * t4xz * (cos(kx) - cos(ky))
    Ham(2, 4) = tvz
    Ham(3, 1) = tvx
    Ham(4, 1) = 2 * t4xz * (cos(kx) - cos(ky))
    Ham(3, 2) = 2 * t4xz * (cos(kx) - cos(ky))
    Ham(4, 2) = tvz
    !---------------------------------------------------------------------
    ! 相互作用矩阵赋值
    U0 = 3.0
    J0 = 0.4
    Umat(1,1) = U0
    Umat(2,2) = U0
    Umat(3,3) = U0
    Umat(4,4) = U0
    Umat(1,2) = J0/2
    Umat(2,1) = J0/2
    Umat(3,4) = J0/2
    Umat(4,3) = J0/2
    !---------------------------------------------------------------------
    ! 单位矩阵
    do k0 = 1,hn
        ones(k0,k0) = 1
    end do
    return
end subroutine
!================================================================================================================================================================================================
subroutine chi0cal(qx,qy,re1)
    ! 计算极化率  返回到re1
    use param
    integer ikx,iky,l1,l2,e1,e2
    real qx,qy,kx,ky
    complex re1(hn,hn)
    do iky = -kn,kn - 1
        ky = pi * iky/kn
        do ikx = -kn,kn - 1
            kx = pi * ikx/kn

            ! k
            call matset(kx,ky)
            call eigSol() 
            valmesh(1,iky,ikx,:) = val(:)
            vecmesh(1,iky,ikx,:,:) = Ham(:,:)

            ! k + q
            call matset(kx + qx,ky + qy)
            call eigSol() 
            valmesh(2,iky,ikx,:) = val(:)
            vecmesh(2,iky,ikx,:,:) = Ham(:,:)
            
            ! 计算极化率
            do l1 = 1,hn   ! orbit ondex
                do l2 = 1,hn
                    do e1 = 1,hn  ! band index
                        do e2 = 1,hn
                            re1(l1,l2) = re1(l1,l2) + (fermi(valmesh(1,iky,ikx,e1)) - fermi(valmesh(2,iky,ikx,e2)))/(im * (omega + 0.0001) + valmesh(1,iky,ikx,e1) - valmesh(2,iky,ikx,e2))&
                                    * conjg(vecmesh(2,iky,ikx,l1,e2)) * vecmesh(2,iky,ikx,l2,e2) * conjg(vecmesh(1,iky,ikx,l2,e1)) * vecmesh(1,iky,ikx,l1,e1)
                        end do
                    end do
                end do
            end do
        end do
    end do
    re1 = re1/(2 * kn)**2
    return 
end subroutine
!================================================================================================================================================================================================
function fermi(ek)
    ! 费米分布函数
    implicit none
    real fermi,ek,kbt
    kbt = 0.001
    fermi = 1/(exp(ek/kbt) + 1)  
    return
end
!================================================================================================================================================================================================
function equivkpq(i0)
    ! 找到BZ中k+q等价与k的索引
    use param
    integer equivkpq,i0
    if (i0 <= kn/2 .and. i0 > -kn/2) equivkpq = i0
    if (i0 > kn/2) equivkpq = i0 - kn
    if (i0 <= -kn/2) equivkpq = i0 + kn
end
!================================================================================================================================================================================================
subroutine eigSol()
    ! 对角化得到本征值w和本征矢量Ham
    use param
    integer m
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','U',hn,Ham,lda,val,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2 * hn + hn**2, int( work( 1 ) ) )
    lrwork = min(1 + 5 * hn + 2 * hn**2, int( rwork( 1 ) ) )
    liwork = min(3 + 5 * hn, iwork( 1 ) )
    call cheevd('V','U',hn,Ham,lda,val,work,lwork,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(11,file = "mes.dat",status = "unknown")
        write(11,*)'The algorithm failed to compute eigenvalues.'
        close(11)
    end if
    ! open(12,file="eigval.dat",status="uknnown")
    ! do m = 1,N
    ! 	write(12,*)m,val(m)
    ! end do
    ! close(12)
    return
end subroutine eigSol
!================================================================================================================================================================================================
subroutine inv(matin,matout)
    ! 矩阵求逆
    use param
    complex,intent(in) :: matin(hn,hn)
    complex:: matout(size(matin,1),size(matin,2))
    real:: work2(size(matin,1))            ! work2 array for LAPACK
    integer::ipiv(size(matin,1))     ! pivot indices
    ! Store matin in matout to prevent it from being overwritten by LAPACK
    matout = matin
    ! SGETRF computes an LU factorization of a general M - by - N matrix A
    ! using partial pivoting with row interchanges .
    call CGETRF(hn,hn,matout,hn,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call CGETRI(hn,matout,hn,ipiv,work2,hn,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
    return
    end subroutine inv
!================================================================================================================================================================================================
```
# MPI并行解读
首先是要在程序中调用MPI的参数库
```shell
use mpi
```
我在代码中的并行思想就是将一个$[-kn,kn-1)$，长度为$2*kn$的循环进行分拆。首先是MPI环境的初始化
```shell
    call MPI_INIT(ierr)     ! 初始化进程
    call MPI_COMM_RANK(MPI_COMM_WORLD, indcore, ierr) ! 得到本进程在通信空间中的rank值,即在组中的逻辑编号(该 rank值为0到p-1间的整数,相当于进程的ID。)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numcore, ierr) !获得进程个数 size, 这里用变量p保存
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    nki = floor(indcore * (2.0 * kn)/numcore) - kn
    nkf = floor((indcore + 1) * (2.0 * kn)/numcore) - kn - 1 
```
{:.success}

并行会同时调用多个CPU进行计算，在调用MPI的时候，每个CPU都会有自己单独的编号，也就是上面的`indcore`这个参数来进行标识，而`numcore`就是在程序执行的时候调用CPU的总数量。此时可以看到对于每个CPU而言，`nki`和`nkf`都是不同的，这样就可以在不同的CPU上面执行$[-kn,kn-1)$这个区间的不同子区间。

这里为了保证并行的效率，需要做到不同的区间之间是没有交叠的，因为后续我们要对不同CPU中计算结果的收集，用的是一个求和的方式，所以如果此时计算区间是有重叠的，那么计算的结果中就会有一部分是重复叠加的，所以一定要保证并行区间是没有重叠，且它们的并集等于$[-kn,kn-1)$.
{:.warning}

下面就是分拆之后的循环在不同的CPU上计算不同的区间
```shell
    do iky = nki,nkf
        qy = pi * iky/kn
        do ikx = -kn,kn - 1
            qx = pi * ikx/kn
            call chi0cal(qx,qy,chi0(iky,ikx,:,:))  ! 得到裸极化率
            call inv(ones - matmul(chi0(iky,ikx,:,:),Umat),temp2)  ! 矩阵求逆
            chi(iky,ikx,:,:) = matmul(temp2,chi0(iky,ikx,:,:))
            temp3 = sum(chi(iky,ikx,:,:))
            local_rechi(iky,ikx,1) = qx
            local_rechi(iky,ikx,2) = qy
            local_rechi(iky,ikx,3) = real(temp3)
            local_rechi(iky,ikx,4) = aimag(temp3)
        end do
    end do
```
每个CPU上都会有自己独立的`local_rechi`变量，下面就是要等所有的CPU计算结束，将结果收集起来
```shell
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_Reduce(local_rechi, rechi, (2 * kn)**2 * 4, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
```
这里调用`call MPI_Barrier(MPI_COMM_WORLD,ierr)`的目的就是等所有的CPU都计算结束，如果有一些CPU已经计算完循环了，而有一些并没有计算结束，这个时候就让程序在这里等候所有的CPU计算都接受，然后调用`call MPI_Reduce(local_rechi, rechi, (2 * kn)**2 * 4, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD,ierr)`来收集数据。其中其中的`rechi`是一个维度和`local_rechi`都完全相同的变量，第三个参数`(2 * kn)**2 * 4`表示要收集的数据量大小。举个例子，你有一个矩阵$A(10,10)$，它一共有100个元素，那么此时需要收集的数据量就是100。第四个参数`MPI_REAL`则用来声明收集数据的类型，还有整型`MPI_INTEGER`和复数型`MPI_COMPLEX`。第五个参数`MPI_SUM`就是数据收集的方式。我这里就需要用求和即可，因为程序在设计的时候不同的CPU计算的结果实际上是存储在矩阵不同位置，没有存储的位置自然就是0，所以用求和就能达到收集数据的目的。第六个参数则是申明在哪个CPU核心上来执行数据收集这个操作，我这里就选择了`root`核心。剩下的就是一些公共参数了，没什么好解释的了。


计算结束之后总是要进行数据存储的，这个操作需要在一个特定的核心上进行
```shell
    if (indcore .eq. 0) then
        char1 = "fortran-chi-"
        write(char2,"(I3.3)")2 * kn
        filename = trim(char1)//trim(char2)
        char1 = ".dat"
        filename = trim(filename)//trim(char1)
        open(12,file = filename)
        do iky = -kn,kn - 1
            do ikx = -kn,kn - 1
                write(12,"(4F8.3)")rechi(iky,ikx,1),rechi(iky,ikx,2),rechi(iky,ikx,3),rechi(iky,ikx,4)
            end do
        end do
        close(12)
    end if
```
程序执行结束之后还申明结束`MPI`并行
```shell
    call MPI_Finalize(ierr)
```

上面实际上只完成了对一个循环的并行，实际情况可以更加丰富。比如我需要用到前面并行得到的结果进行后续计算，而且后续的计算也需要并行，此时就需要对计算得到的变量进行数据广播，然后再进行并行计算，这个方法在后续会继续更新，目前的代码还没有涉及到这个问题。


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