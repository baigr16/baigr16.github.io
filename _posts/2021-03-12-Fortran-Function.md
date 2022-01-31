---
title: Fortran 常用函数及操作
tags: Code Fortran 
layout: article
license: true
toc: true
key: a20210312
pageview: true
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
因为自己在学习中最常使用的是Fortran，有时候要用到mkl库中的一些函数，但是这些函数的调用参数有很多，所以将自己常用的一些记录下来，以后忘了可以
快速的查阅。
<!--more-->
# 1. cheevd复厄密矩阵对角化

> 厄密矩阵本征值与本征矢的求解

# 2.getrf,getri 矩阵求逆

> getrf 对一个矩阵进行LU分解

> getri 计算由LU分解后矩阵的逆 

```fortran
program main
use lapack95	!调用函数库，确保使用的函数存在
!如果时用的是intel fortran，则在编译的时候注释上一行，加入-mkl编译选项即可。
implicit none
integer, parameter :: n = 3
integer :: i, j, temp(n)
real(kind=8) :: a(n,n), aa(n,n)

call random_seed()
call random_number(a)
aa = a

write(*,'(1x,a)') "a = "
do i = 1, n
        write(*,'(*(f12.6,3x))') a(i,:)
end do

!使用库函数求逆
call getrf( a, temp )!先进性LU分解
call getri( a, temp )!然后进行矩阵求逆

write(*,'(1x,a)') 'inv(a) = '
do i = 1, n
        write(*,'(*(f12.6,3x))') a(i,:)
end do

write(*,'(1x,a)') "checking..."
aa = matmul(aa,a)  !// 原矩阵与其逆矩阵的结果为单位矩阵
do i = 1, n
        write(*,'(*(f12.6,3x))') aa(i,:)
end do

end program main
```

结果如下

```fortran
a =
    0.959299       0.268247       0.274620
    0.013673       0.082084       0.275984
    0.056097       0.730892       0.177709
 inv(a) =
    1.072178      -0.876915      -0.295021
   -0.074784      -0.888503       1.495423
   -0.030876       3.931107      -0.430155
 checking...
    1.000000       0.000000      -0.000000
   -0.000000       1.000000       0.000000
   -0.000000       0.000000       1.000000
```



# 3.?geev一般矩阵对角化
```fortran
call sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
call dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
call cgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)  
```

# 矩阵直积
最近遇到一些要批量计算的问题,其中涉及到许多Pauli矩阵的直积,这里就整理一下如何利用Fortran来计算矩阵直积
```fortran
! Author:YuXuanLi
! E-Mail:yxli406@gmail.com
    module pub
    implicit none
    integer xn,yn,len2,N
    parameter(xn = 30,yn = 30,N = xn*yn*4,len2 = xn*yn)
    complex,parameter::im = (0.0,1.0) 
    real,parameter::pi = 3.14159265359
    complex Ham(N,N) 
    integer bry(8,len2)
    complex sx(2,2),sy(2,2),sz(2,2),s0(2,2)k   
    integer info
    end module pub
!========== PROGRAM START ==========================
    program sol
    use pub
    call test()
    stop
    end program sol
!===================================================
    Module Kron_Mod
    Implicit None
    contains
        Subroutine Kron( A , B , H )
        complex, Intent( IN )  :: A(:,:) , B(:,:)
        complex, Intent( OUT ) :: H(:,:)
        Integer :: i , j , m , n , p , q
        m = size( A , dim = 1 )
        n = size( A , dim = 2 )
        p = size( B , dim = 1 )
        q = size( B , dim = 2 )
        Do i = 1 , m
            Do j = 1 , n
            H( p*(i-1)+1 : p*i , q*(j-1)+1 : q*j ) = B * A(i,j)
            End Do
        End Do
        End Subroutine Kron
    
    End Module Kron_Mod
!=========================================================
    ! Build module in here, with Block Matrix operations. 
    Module Kronecker
        contains
  ! Takes in Matrices A(i,j),B(k,l), assumed 2D, returns Kronecker Product C(i*k,j*l)
        function KronProd(A,B) result(C)
         IMPLICIT NONE
         real, dimension (:,:), intent(in)  :: A, B
         real, dimension (:,:), allocatable :: C
         integer :: i = 0, j = 0, k = 0, l = 0
         integer :: m = 0, n = 0, p = 0, q = 0
         allocate(C(size(A,1)*size(B,1),size(A,2)*size(B,2)))
         C = 0
         do i = 1,size(A,1)
            do j = 1,size(A,2)
                n=(i-1)*size(B,1) + 1
                m=n+size(B,1) - 1
                p=(j-1)*size(B,2) + 1
                q=p+size(B,2) - 1
                C(n:m,p:q) = A(i,j)*B
            end do
         end do
        end function KronProd  
    !----------------------------------------------------------------  
  ! Takes in Matrices A(i,j),B(k,l), assumed 2D, returns Direct sum
  ! C(i+k,j+l)
        function DirSum(A,B) result(C)
         real, dimension (:,:), intent(in)  :: A, B
         real, dimension (:,:), allocatable :: C
         integer :: p = 0, q = 0
         allocate(C(size(A,1)+size(B,1),size(A,2)+size(B,2)))
         C = 0
         p = size(A,1) + size(B,1) 
         q = size(A,2) + size(B,2) 
         C(1:size(A,1),1:size(A,2)) = A
         C(size(A,1)+1:p,size(A,2)+1:q) = B
         return
         end function DirSum
    !--------------------------------------------------------
  ! Takes 2 vectors, A(i),B(j), returns Direct Sum C(i+j)
        function VecDirSum(A,B) result(C)
         real, dimension (:), intent(in)  :: A, B
         real, dimension (:), allocatable :: C
         allocate(C(size(A)+size(B)))
         C = 0
         C(1:size(A)) = A
         C(size(A)+1:size(A)+size(B)) = B
         return
         end function VecDirSUm
    end module Kronecker
!===================================================
    subroutine test()
    use pub
    use Kron_Mod
    complex H1(2**2,2**2),H2(2**2,2**2)
    call Pauli()
    call Kron(sx,sx,H)
    H2 = KronProd(sx,sx)
    Do i = 1 , 2**2
        Write(*,'(4(f5.1,1x))') H1( :, i)
    end do
    end subroutine test
!===========================================================
    subroutine Pauli()
    use pub
    sx(1,2) = 1
    sx(2,1) = 1
    !----
    sy(1,2) = -im
    sy(2,1) = im
    !-----
    sz(1,1) = 1
    sz(2,2) = -1
    !-----
    s0(1,1) = 1
    s0(2,2) = 1
    end subroutine Pauli
```

这里提供了两个模块,里面都有计算矩阵直积的函数,但有点小问题,那就是只能进行实数矩阵的直积,复数矩阵的直积暂时还没有解决.
{:.warning}

这里时网上的另外一种计算程序
```fortran
Module Kron_Mod
Implicit None
Integer , parameter , private :: DP = Selected_Real_Kind( 9 )
contains

    Subroutine Kron( A , B , H )
    Real( Kind = DP ) , Intent( IN )  :: A(:,:) , B(:,:)
    Real( Kind = DP ) , Intent( OUT ) :: H(:,:)
    Integer :: i , j , m , n , p , q
    m = size( A , dim = 1 )
    n = size( A , dim = 2 )
    p = size( B , dim = 1 )
    q = size( B , dim = 2 )
    Do i = 1 , m
        Do j = 1 , n
        H( p*(i-1)+1 : p*i , q*(j-1)+1 : q*j ) = B * A(i,j)
        End Do
    End Do
    End Subroutine Kron

End Module Kron_Mod
!========================================================
Program www_fcode_cn
use Kron_Mod
Implicit None
Integer , parameter :: DP = Selected_Real_Kind( 9 )
Integer , parameter :: m=2 , n=3 , p=4, q=5 , index1 = m*p , index2 = n*q
Real(kind=DP) :: A(m,n) , B(p,q) , H(index1,index2),sx(m,m),H2(m**2,m**2)
integer :: i
sx(1,2) = 1
sx(2,1) = 1
A = reshape( (/1,2,3,4,5,6/) , (/2,3/) )
B = reshape( (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20/) , (/4,5/) )
call Kron( A , B , H )
call Kron( sx , sx , H2 )
Do i = 1 , m**2
! Write(*,'(8(f5.1,1x))') H( :, i)
Write(*,'(4(f5.1,1x))') H2( :, i)
End do
End Program www_fcode_cn
```
# 矩阵求逆
```fortran
subroutine inv(ndim,Amat)
    ! ndim is dimension of Amat,the input Amat will be written which is result
    implicit none
    integer,parameter::dp = 8
    integer::i
    integer::info
    integer,intent(in)::ndim   ! in 代表这个值只能是个输入值,不会被改变
!    IPIV   : INTEGER. Array,DIMENSION at least max(1,n). The pivot indices that define
!    the permutation matrix P; row i of the matrix was interchanged with
!    row ipiv(i). Corresponds to the single precision factorization (if info=
!    0 and iter ≥ 0) or the double precision factorization (if info= 0 and
!    iter < 0).
    integer,allocatable::ipiv(:)
    complex(dp),parameter::zone = (1.0d0,0.0d0)
!    Amat  :
!    Overwritten by the factors L and U from the factorization of A = P*L*U;
!    the unit diagonal elements of L are not stored.
!    If iterative refinement has been successfully used (info= 0 and iter≥
!    0),then A is unchanged.
!    If double precision factorization has been used (info= 0 and iter <
!    0),then the array A contains the factors L and U from the factorization
!    A = P*L*U; the unit diagonal elements of L are not stored.
    complex(dp),intent(inout):: Amat(ndim,ndim)
!    Bmat  :
!    Overwritten by the solution matrix X for dgesv,sgesv,zgesv,zgesv.
    complex(dp),allocatable::Bmat(:,:)
    allocate(ipiv(ndim))
    allocate(Bmat(ndim,ndim))
    ipiv=0
    ! unit matrix
    Bmat= (0d0,0d0)
    do i=1,ndim
        Bmat(i,i)= zone
    end do
    call zgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)
    if(info.ne.0)print *,'something wrong with zgesv'
    Amat=Bmat
    return
    end subroutine inv
```

# 矩阵相乘
```fortran
! performs matrix-matrix multiply
! C := alpha*op( A )*op( B ) + beta*C
    subroutine mat_mul(nmatdim,A,B,C)  
    ! nmatdim is dimension of matrix
    implicit none
    integer,parameter::dp = kind(1.0d0) ! double precision(精度控制)
    integer,intent(in)::nmatdim    
    complex(dp)::ALPHA
    complex(dp)::BETA 
    complex(dp),intent(in) ::A(nmatdim ,nmatdim)
    complex(dp),intent(in) ::B(nmatdim ,nmatdim)
    !complex(dp)::mat_mul(nmatdim,nmatdim)
    complex(dp),intent(out)::C(nmatdim,nmatdim)
    alpha = 1.0d0 
    beta = 0.0d0
    C(:,:) = (0.0d0,0.0d0)
    call ZGEMM('N','N',nmatdim,nmatdim,nmatdim,ALPHA,A,nmatdim,B,nmatdim,BETA,C,nmatdim)
    return
    end subroutine mat_mul
```

# 读取文件行数
通常在读取文件的时候, 并不会指导文件到底有多少行, 这个子过程就是用来确定在读取到文件结尾的时候终止循环, 从而指导一共有多少行数据
```fortran
subroutine main2()
! 读取不明行数的文件
implicit none
integer count,stat
real h1,h2,h3
open(1,file = "wavenorm.dat")
do while (.true.)
    count = count + 1
    read(1,*,iostat = STAT)h1,h2,h3
    if(stat .ne. 0) exit ! 当这个参数不为零的时候,证明读取到文件结尾
end do
write(*,*)h1,h2,h3
write(*,*)count
close(1)
return
end subroutine 
```

# 读取不明行数的文件
```fortran
    program main
    implicit none
    integer m1,m2,m3
    call main2()
    stop
    end program 
!=======================================================    
    subroutine main2()
    ! 读取不明行数的文件
    implicit none
    integer count,stat
    real h1,h2,h3,h22
    h1 = 0
    h2 = 0
    h3 = 0
    h22 = 0
    open(1,file = "da.dat")
    open(2,file = "da2.dat")
    count = 0
    do while (.true.)
        count = count + 1
        h22 = h2
        read(1,*,iostat = STAT)h1,h2,h3
        if(h22.ne.h2)write(2,*)""  ! 在这里加空行是为了gnuplot绘制密度图
        write(2,999)h1,h2,h3   ! 数据格式化
        if(stat .ne. 0) exit ! 当这个参数不为零的时候,证明读取到文件结尾
    end do
    write(*,*)h1,h2,h3
    write(*,*)count
    close(1)
    close(2)
999 format(3f11.6)
    return
    end subroutine 
```

# 整数/浮点数转换为字符串
```fortran
    subroutine plot(m3)
    ! 将能带图绘制到一起
    use pub
    character*20::str1,str2,str3,str4,str5,str6
    integer m3
    str1 = "half-ud"
    write(str2,"(I2.2)")m3 将整数转换成字符串
    str3 = ".gnu"
    str4 = trim(str1)//trim(str2)//trim(str3)  ! trim()在拼接字符串的时候删除字符串中多余的空格
    write(str5,'(f5.2)')dx ! 将一个小数转成字符串
    write(*,*)str2
    write(*,*)str5
    return
    end subroutine plot
```