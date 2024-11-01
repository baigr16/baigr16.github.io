---
title: Fortran常用对角化函数的封装
tags:  Fortran Code 
layout: article
license: true
toc: true
key: a20240420
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
这里整理一下在用Fortran编写物理程序的时候经常用的几个函数封装，目前只是简单的几个矩阵对角化函数，后续有新的函数也会继续增加。
{:.info}
<!--more-->
# 前言
最近因为遇到了计算瓶颈，需要用`Fortran`来重写代码，其中就遇到了对矩阵的一些操作，虽然[lapack](https://www.netlib.org/lapack/)已经提供了很好的函数，但是在实际程序编写中还是有些麻烦，这里就对目前我用到的函数进行了一些封装，方便调用。

## 一般方矩阵对角化
```fortran
subroutine diagonalize_general_matrix(matdim,matin,matout,mateigval)
    use param,only:dp,im
    implicit none
    integer matdim, LDA, LDVL, LDVR, LWMAX, INFO, LWORK
    complex(dp) mateigval(matdim),matin(matdim,matdim),matout(matdim,matdim)
    real,allocatable::VL(:,:),VR(:,:),WR(:),WI(:),WORK(:)
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    allocate(VL( LDVL, matdim), VR(LDVL,matdim),  WR(matdim), WI(matdim), WORK(LWMAX))
    LWORK = -1
    matout = matin
    CALL SGEEV( 'V', 'V', matdim, matout, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    CALL SGEEV( 'V', 'V', matdim, matout, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The algorithm failed to compute eigenvalues.'
        STOP
     END IF
     do info = 1,matdim
        mateigval(info) = WR(info) + im * WI(info)
    end do
    matout = VL  ! 将左矢量输出
    return
end subroutine
```
## 一般方矩阵对角化(双精度)
```fortran
subroutine diagonalize_general_matrix(matdim, matin, matout, mateigval)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, intent(in) :: matdim
    integer :: LDA, LDVL, LDVR, LWMAX, INFO, LWORK
    complex(dp), intent(out) :: mateigval(matdim), matout(matdim, matdim)
    complex(dp), intent(in) :: matin(matdim, matdim)
    real(dp), allocatable :: VL(:,:), VR(:,:), WR(:), WI(:), WORK(:)
    ! 矩阵维度参数设置
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    matout = matin
    ! 分配工作空间
    allocate(VL(LDVL, matdim), VR(LDVR, matdim), WR(matdim), WI(matdim), WORK(LWMAX))
    ! 初次调用 DGEEV 获取最佳工作空间大小
    LWORK = -1
    CALL DGEEV('V', 'V', matdim, matout, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    ! 根据第一次调用的返回值重新调整工作空间大小
    LWORK = MIN(LWMAX, INT(WORK(1)))
    deallocate(WORK)
    allocate(WORK(LWORK))
    ! 第二次调用 DGEEV 计算特征值和特征向量
    CALL DGEEV('V', 'V', matdim, matout, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    ! 检查对角化是否成功
    IF (INFO > 0) THEN
        WRITE(*,*) 'The algorithm failed to compute eigenvalues.'
        STOP
    END IF
    ! 将实部和虚部合成为复数本征值
    do info = 1, matdim
        mateigval(info) = WR(info) + (0.0, 1.0) * WI(info)
    end do
    ! 将右特征向量存储在 matout 中
    matout = VR
    ! 释放工作空间
    deallocate(VL, VR, WR, WI, WORK)
    return
end subroutine diagonalize_general_matrix
```

## 复数矩阵对角化(单精度)
```fortran
subroutine diagonalize_complex_matrix(matdim,matin,matout,mateigval)
    ! 对角化一般复数矩阵,这里的本征值是个复数
    ! matin 输入矩阵   matout 本征矢量    mateigval  本征值
    integer matdim,LDA,LDVL,LDVR,LWMAX,INFO,LWORK
    complex,intent(in)::matin(matdim,matdim)
    complex,intent(out)::matout(matdim,matdim)
    complex,intent(out)::mateigval(matdim)
    REAL,allocatable::RWORK(:)
    complex,allocatable::WORK(:)
    complex,allocatable::VL(:,:)
    complex,allocatable::VR(:,:)
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    ! write(*,*)matin
    allocate(RWORK(2 * matdim))
    allocate(VL(LDVL,matdim))
    allocate(VR(LDVR, matdim))
    allocate(WORK(LWMAX))
    matout = matin

    LWORK = -1
    call cgeev( 'V', 'N', matdim, matout, LDA, mateigval, VL, LDVL,VR, LDVR, WORK, LWORK, RWORK, INFO)
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    call cgeev( 'V', 'N', matdim, matout, LDA, mateigval, VL, LDVL,VR, LDVR, WORK, LWORK, RWORK, INFO )
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    !    STOP
    END IF
    matout = VL  ! 将左矢量输出
    return
end subroutine diagonalize_complex_matrix
```

## 复数矩阵对角化(双精度)
```fortran
subroutine diagonalize_complex_matrix(matdim, matin, matout, mateigval)
    ! 对角化一般复数矩阵，这里的本征值是复数(双精度)
    ! matin 输入矩阵, matout 本征矢量, mateigval 本征值
    implicit none
    integer, intent(in) :: matdim
    integer :: LDA, LDVL, LDVR, LWMAX, INFO, LWORK
    integer, parameter :: dp = kind(1.0d0)
    complex(dp), intent(in) :: matin(matdim, matdim)
    complex(dp), intent(out) :: matout(matdim, matdim)
    complex(dp), intent(out) :: mateigval(matdim)
    real(dp), allocatable :: RWORK(:)
    complex(dp), allocatable :: WORK(:)
    complex(dp), allocatable :: VL(:,:)
    complex(dp), allocatable :: VR(:,:)
    ! 设置矩阵维度参数
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    ! 分配工作空间
    allocate(RWORK(2 * matdim))
    allocate(VL(LDVL, matdim))
    allocate(VR(LDVR, matdim))
    allocate(WORK(LWMAX))
    ! 将输入矩阵赋值给 matout
    matout = matin
    ! 初次调用 ZGEEV 获取最佳工作空间大小
    LWORK = -1
    call ZGEEV('V', 'V', matdim, matout, LDA, mateigval, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    ! 根据第一次调用的返回值重新调整工作空间大小
    LWORK = min(LWMAX, int(WORK(1)))
    deallocate(WORK)
    allocate(WORK(LWORK))
    ! 第二次调用 ZGEEV 计算本征值和本征矢量
    call ZGEEV('V', 'V', matdim, matout, LDA, mateigval, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    ! 检查对角化是否成功
    if (INFO > 0) then
        write(*, *) 'The algorithm failed to compute eigenvalues.'
    end if
    ! 将右特征向量存储在 matout 中
    matout = VR
    ! 释放工作空间
    deallocate(WORK)
    deallocate(RWORK)
    deallocate(VL)
    deallocate(VR)
    return
end subroutine diagonalize_complex_matrix
```


## 矩阵求逆(单精度)
```fortran
subroutine Matrix_Inv(matdim,matin,matout)
    ! 矩阵求逆 
    implicit none
    integer matdim,dp,info
    parameter(dp = kind(1.0))
    complex(dp),intent(in) :: matin(matdim,matdim)
    complex(dp):: matout(size(matin,1),size(matin,2))
    real(dp):: work2(size(matin,1))            ! work2 array for LAPACK
    integer::ipiv(size(matin,1))     ! pivot indices
    ! Store matin in matout to prevent it from being overwritten by LAPACK
    matout = matin
    ! SGETRF computes an LU factorization of a general M - by - N matrix A
    ! using partial pivoting with row interchanges .
    call CGETRF(matdim,matdim,matout,matdim,ipiv,info)
    ! if (info.ne.0) stop 'Matrix is numerically singular!'
    if (info.ne.0)  write(*,*)'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call CGETRI(matdim,matout,matdim,ipiv,work2,matdim,info)
    ! if (info.ne.0) stop 'Matrix inversion failed!'
    if (info.ne.0) write(*,*)'Matrix inversion failed!'
    return
end subroutine Matrix_Inv
```
## 矩阵求逆(双精度)
```fortran
subroutine Matrix_Inv(matdim, matin, matout)
    ! 矩阵求逆(双精度)
    implicit none
    integer, intent(in) :: matdim
    integer :: info
    integer, parameter :: dp = kind(1.0d0)
    complex(dp), intent(in) :: matin(matdim, matdim)
    complex(dp) :: matout(matdim, matdim)
    complex(dp), allocatable :: work2(:)       ! 工作空间数组
    integer :: ipiv(matdim)                     ! pivot indices
    ! 将输入矩阵 matin 复制到 matout，防止 LAPACK 函数覆盖 matin
    matout = matin
    ! 确定工作空间大小
    allocate(work2(matdim))
    ! ZGETRF 用于矩阵 LU 分解
    call ZGETRF(matdim, matdim, matout, matdim, ipiv, info)
    if (info /= 0) then
        write(*, *) 'Matrix is numerically singular!'
        return
    end if
    ! ZGETRI 用于求解逆矩阵
    call ZGETRI(matdim, matout, matdim, ipiv, work2, matdim, info)
    if (info /= 0) then
        write(*, *) 'Matrix inversion failed!'
        return
    end if
    ! 释放工作空间
    deallocate(work2)
    return
end subroutine Matrix_Inv
```


## 厄米矩阵对角化(单精度)
```fortran
subroutine diagonalize_Hermitian_matrix(matdim,matin,matout,mateigval)
    !  厄米矩阵对角化
    ! matin 输入矩阵   matout 本征矢量    mateigval  本征值
    integer matdim
    integer lda0,lwmax0,lwork,lrwork,liwork,info
    complex matin(matdim,matdim),matout(matdim,matdim)
    real mateigval(matdim)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    !-----------------
    lda0 = matdim
    lwmax0 = 2 * matdim + matdim**2
    allocate(work(lwmax0))
    allocate(rwork(1 + 5 * matdim + 2 * matdim**2))
    allocate(iwork(3 + 5 * matdim))
    matout = matin
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','U',matdim,matout,lda0,mateigval,work,lwork ,rwork,lrwork,iwork,liwork,info)
    lwork = min(2 * matdim + matdim**2, int( work( 1 ) ) )
    lrwork = min(1 + 5 * matdim + 2 * matdim**2, int( rwork( 1 ) ) )
    liwork = min(3 + 5 * matdim, iwork( 1 ) )
    call cheevd('V','U',matdim,matout,lda0,mateigval,work,lwork,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(11,file = "mes.dat",status = "unknown")
        write(11,*)'The algorithm failed to compute eigenvalues.'
        close(11)
    end if
    return
end subroutine diagonalize_Hermitian_matrix
```

## 厄米矩阵对角化(双精度)
```fortran
subroutine diagonalize_Hermitian_matrix(matdim, matin, matout, mateigval)
    ! 厄米矩阵对角化
    ! matin 输入矩阵, matout 本征矢量, mateigval 本征值
    implicit none
    integer, parameter :: dp = kind(1.0d0)  ! 双精度
    integer, intent(in) :: matdim
    integer :: lda0, lwmax0, lwork, lrwork, liwork, info
    complex(dp), intent(in) :: matin(matdim, matdim)
    complex(dp), intent(out) :: matout(matdim, matdim)
    real(dp), intent(out) :: mateigval(matdim)
    complex(dp), allocatable :: work(:)
    real(dp), allocatable :: rwork(:)
    integer, allocatable :: iwork(:)

    !-----------------
    lda0 = matdim
    lwmax0 = 2 * matdim + matdim**2
    allocate(work(lwmax0))
    allocate(rwork(1 + 5 * matdim + 2 * matdim**2))
    allocate(iwork(3 + 5 * matdim))

    matout = matin
    lwork = -1
    liwork = -1
    lrwork = -1

    ! 初次调用 zheevd 获取最佳工作空间大小
    call zheevd('V', 'U', matdim, matout, lda0, mateigval, work, lwork, rwork, lrwork, iwork, liwork, info)

    ! 根据第一次调用的返回值重新调整工作空间大小
    lwork = min(2 * matdim + matdim**2, int(work(1)))
    lrwork = min(1 + 5 * matdim + 2 * matdim**2, int(rwork(1)))
    liwork = min(3 + 5 * matdim, iwork(1))

    ! 重新分配工作空间
    deallocate(work)
    allocate(work(lwork))
    deallocate(rwork)
    allocate(rwork(lrwork))
    deallocate(iwork)
    allocate(iwork(liwork))

    ! 第二次调用 zheevd 计算本征值和本征矢量
    call zheevd('V', 'U', matdim, matout, lda0, mateigval, work, lwork, rwork, lrwork, iwork, liwork, info)

    ! 错误处理
    if (info > 0) then
        open(11, file = "mes.dat", status = "unknown")
        write(11, *) 'The algorithm failed to compute eigenvalues.'
        close(11)
    end if

    return
end subroutine diagonalize_Hermitian_matrix
```
**提示**

如果使用了双精度，那么在MPI并行的时候，数据收集也需要对应的修改为双精度
{:.info}
```shell
    call MPI_Barrier(MPI_COMM_WORLD,ierr)   ! 等所有核心都计算完成
    call MPI_Reduce(U0_mpi, U0_list, Un, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
    call MPI_Reduce(da_mpi, da_list, Un, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
    call MPI_Reduce(Ds_con_mpi, Ds_con_list, 4 * Un, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
```

## 实数矩阵对角化(单精度)
```fortran
subroutine diagonalize_real_matrix(matdim,matin,matout,mateigval)
    ! 对角化一般复数矩阵,这里的本征值是个复数
    ! matin 输入矩阵   matout 本征矢量    mateigval  本征值
    integer matdim,LDA,LDVL,LDVR,LWMAX,INFO,LWORK
    real,intent(in)::matin(matdim,matdim)
    real,intent(out)::matout(matdim,matdim)
    complex,intent(out)::mateigval(matdim)
    real valre(matdim),valim(matdim)
    real,allocatable::WORK(:)
    real,allocatable::VL(:,:)
    real,allocatable::VR(:,:)
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    ! write(*,*)matin
    allocate(VL(LDVL,matdim))
    allocate(VR(LDVR, matdim))
    allocate(WORK(LWMAX))
    matout = matin

    LWORK = -1
    call sgeev( 'V', 'N', matdim, matout, LDA, valre,valim, VL, LDVL,VR, LDVR, WORK, LWORK, INFO)
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    call sgeev( 'V', 'N', matdim, matout, LDA, valre,valim, VL, LDVL,VR, LDVR, WORK, LWORK, INFO)
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    !    STOP
    END IF
    do info = 1,matdim
        mateigval(info) = valre(info) + im * valim(info)
    end do
    matout = VL  ! 将左矢量输出
    return
end subroutine diagonalize_real_matrix
```
实数矩阵的双精度版本实际上就是一般矩阵对角化

# 完整代码&测试
```fortran
module param
    implicit none
    integer, parameter :: dp = kind(1.0)
    integer num_wann,nrpts,hn
    parameter(hn = 4)
    real(dp) ones(hn,hn)
    real(dp),parameter::pi = 3.1415926535897
    complex(dp),parameter::im = (0.,1.) !Imagine unit
    real(dp) tsig,tpi,tsig2,tpi2,mu  ! 哈密顿量参数
    parameter(tsig = 2.0,tpi = tsig/1.56,tsig2 = tsig/7.0,tpi2 = tsig2/1.56)
    complex(dp),allocatable::HmnR(:,:,:)
    real(dp),allocatable::ndegen(:)
    real(dp),allocatable::irvec(:,:)
    integer,allocatable::BZklist(:,:)
end module
!==============================================================================================================================================================
program main
    use param
    implicit none
    complex(dp) Ham(hn,hn),temp1(hn,hn),temp2(hn),A1(4,4)
    real(dp) temp3(hn),A0(10,10)
    character(len = 20) char1,char2,char3
    integer i1,i2

     do i1 = 1,10
         do i2 = 1,10
             A0(i1,i2) = i1 + i2
         end do
     end do
     char1 = "("
     write(char2,"(I4)")10
     char3 = "F8.4)"
     char2 = trim(char1)//trim(char2)//trim(char3)
    !  write(*,*)char2
     write(*,"(10F8.4)")A0
     write(*,*)"--------------------------------------"
     write(*,char2)A0
    
     call matset(0.3,0.2,Ham)
     write(*,*)Ham
     write(*,*)"-----------------------------------------"
     call diagonalize_complex_matrix(4,Ham,temp1,temp2)
     call diagonalize_Hermitian_matrix(4,Ham,temp1,temp3)
     write(*,*)temp2
     write(*,*)temp3
     call Matrix_Inv(4,Ham,temp1)
     write(*,*)matmul(Ham,temp1)

     call diagonalize_real_matrix(4,A0,A1,temp2)
     write(*,*)temp2

    stop
end program main
!==============================================================================================================================================================
subroutine matset(kx,ky,Ham)
    ! 矩阵赋值,返回Ham
    use param
    implicit none
    real(dp) kx,ky
    integer k0
    complex(dp) Ham(hn,hn)
    complex(dp) h1(hn,hn),h2(hn,hn),h13,h14,h23,h24,hn11,hn22,hn12
    
    h13 = 1.0/4*(tpi + 3.0 * tsig)*(exp(im*(ky/(2 * sqrt(3.0))- kx/2.0)) + exp(im*(ky/(2 * sqrt(3.0)) + kx/2.0)) ) + tpi*exp(-im * ky/sqrt(3.0))
    h14 = -sqrt(3.0)/4 * (tpi - tsig) * (-1 + exp(im * kx)) * exp(-1.0/6.0 * im * (3.0 * kx - sqrt(3.0) * ky ))
    h23 = h14
    h24 = 1.0/4 * (3.0 * tpi + tsig) * ( exp(im * (ky/(2 * sqrt(3.0)) - kx/2.0 ) ) +  exp(im * (ky/(2 * sqrt(3.0)) + kx/2.0 ) )) + tsig * exp(-im * ky/sqrt(3.0))
    hn11 = (3.0 * tpi2 + tsig2) * cos(kx/2.0) * cos(sqrt(3.0) * ky/2.0) + 2 * tsig2 * cos(kx)
    hn12 = sqrt(3.0) *   (tpi2 - tsig2) * sin(kx/2.0) * sin(sqrt(3.0) * ky /2.0)
    hn22 = (tpi2 + 3.0 * tsig2) * cos(kx/2.0) * cos(sqrt(3.0) * ky / 2) + 2 * tpi2 * cos(kx)

    ! 单位矩阵
    do k0 = 1,hn
        ones(k0,k0) = 1
    end do

    
    
    H1(1,3) = h13
    H1(3,1) = conjg( H1(1,3))
    H1(1,4) = h14
    H1(4,1) = conjg(H1(1,4))
    H1(2,3) = h23
    H1(3,2) = conjg(H1(2,3))
    H1(2,4) = h24
    H1(4,2) = conjg(H1(2,4))
    !-------------
    H2(1,1) = hn11
    H2(1,2) = hn12
    H2(2,1) = conjg(H2(1,2))
    H2(2,2) = hn22
    H2(3,3) = hn11
    H2(3,4) = hn12
    H2(4,3) = conjg(H2(3,4))
    H2(4,4) = hn22
    !-------------
    Ham = 0.0
    Ham = H1 + H2 - mu * ones
    !---------------------------------------------------------------------
    return
end subroutine

!==============================================================================================================================================================
subroutine diagonalize_complex_matrix(matdim,matin,matout,mateigval)
    ! 对角化一般复数矩阵,这里的本征值是个复数
    ! matin 输入矩阵   matout 本征矢量    mateigval  本征值
    integer matdim,LDA,LDVL,LDVR,LWMAX,INFO,LWORK
    complex,intent(in)::matin(matdim,matdim)
    complex,intent(out)::matout(matdim,matdim)
    complex,intent(out)::mateigval(matdim)
    REAL,allocatable::RWORK(:)
    complex,allocatable::WORK(:)
    complex,allocatable::VL(:,:)
    complex,allocatable::VR(:,:)
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    ! write(*,*)matin
    allocate(RWORK(2 * matdim))
    allocate(VL(LDVL,matdim))
    allocate(VR(LDVR, matdim))
    allocate(WORK(LWMAX))
    matout = matin

    LWORK = -1
    call cgeev( 'V', 'N', matdim, matout, LDA, mateigval, VL, LDVL,VR, LDVR, WORK, LWORK, RWORK, INFO)
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    call cgeev( 'V', 'N', matdim, matout, LDA, mateigval, VL, LDVL,VR, LDVR, WORK, LWORK, RWORK, INFO )
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    !    STOP
    END IF
    matout = VL  ! 将左矢量输出
    return
end subroutine diagonalize_complex_matrix
!================================================================================================================================================================================================
subroutine Matrix_Inv(matdim,matin,matout)
    ! 矩阵求逆 
    implicit none
    integer matdim,dp,info
    parameter(dp = kind(1.0))
    complex(dp),intent(in) :: matin(matdim,matdim)
    complex(dp):: matout(size(matin,1),size(matin,2))
    real(dp):: work2(size(matin,1))            ! work2 array for LAPACK
    integer::ipiv(size(matin,1))     ! pivot indices
    ! Store matin in matout to prevent it from being overwritten by LAPACK
    matout = matin
    ! SGETRF computes an LU factorization of a general M - by - N matrix A
    ! using partial pivoting with row interchanges .
    call CGETRF(matdim,matdim,matout,matdim,ipiv,info)
    ! if (info.ne.0) stop 'Matrix is numerically singular!'
    if (info.ne.0)  write(*,*)'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call CGETRI(matdim,matout,matdim,ipiv,work2,matdim,info)
    ! if (info.ne.0) stop 'Matrix inversion failed!'
    if (info.ne.0) write(*,*)'Matrix inversion failed!'
    return
end subroutine Matrix_Inv
!================================================================================================================================================================================================
subroutine diagonalize_Hermitian_matrix(matdim,matin,matout,mateigval)
    !  厄米矩阵对角化
    ! matin 输入矩阵   matout 本征矢量    mateigval  本征值
    integer matdim
    integer lda0,lwmax0,lwork,lrwork,liwork,info
    complex matin(matdim,matdim),matout(matdim,matdim)
    real mateigval(matdim)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    !-----------------
    lda0 = matdim
    lwmax0 = 2 * matdim + matdim**2
    allocate(work(lwmax0))
    allocate(rwork(1 + 5 * matdim + 2 * matdim**2))
    allocate(iwork(3 + 5 * matdim))
    matout = matin
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','U',matdim,matout,lda0,mateigval,work,lwork ,rwork,lrwork,iwork,liwork,info)
    lwork = min(2 * matdim + matdim**2, int( work( 1 ) ) )
    lrwork = min(1 + 5 * matdim + 2 * matdim**2, int( rwork( 1 ) ) )
    liwork = min(3 + 5 * matdim, iwork( 1 ) )
    call cheevd('V','U',matdim,matout,lda0,mateigval,work,lwork,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(11,file = "mes.dat",status = "unknown")
        write(11,*)'The algorithm failed to compute eigenvalues.'
        close(11)
    end if
    return
end subroutine diagonalize_Hermitian_matrix
!==============================================================================================================================================================
subroutine diagonalize_real_matrix(matdim,matin,matout,mateigval)
    ! 对角化一般复数矩阵,这里的本征值是个复数
    ! matin 输入矩阵   matout 本征矢量    mateigval  本征值
    integer matdim,LDA,LDVL,LDVR,LWMAX,INFO,LWORK
    real,intent(in)::matin(matdim,matdim)
    real,intent(out)::matout(matdim,matdim)
    complex,intent(out)::mateigval(matdim)
    real valre(matdim),valim(matdim)
    real,allocatable::WORK(:)
    real,allocatable::VL(:,:)
    real,allocatable::VR(:,:)
    LDA = matdim
    LDVL = matdim
    LDVR = matdim
    LWMAX = 2 * matdim + matdim**2
    ! write(*,*)matin
    allocate(VL(LDVL,matdim))
    allocate(VR(LDVR, matdim))
    allocate(WORK(LWMAX))
    matout = matin

    LWORK = -1
    call sgeev( 'V', 'N', matdim, matout, LDA, valre,valim, VL, LDVL,VR, LDVR, WORK, LWORK, INFO)
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    call sgeev( 'V', 'N', matdim, matout, LDA, valre,valim, VL, LDVL,VR, LDVR, WORK, LWORK, INFO)
    IF( INFO.GT.0 ) THEN
       WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    !    STOP
    END IF
    do info = 1,matdim
        mateigval(info) = valre(info) + im * valim(info)
    end do
end subroutine diagonalize_real_matrix
```

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