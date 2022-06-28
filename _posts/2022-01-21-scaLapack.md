---
title: Fortran稀疏矩阵并行化求解
tags: Study Fortran
layout: article
license: true
toc: true
key: a20220121
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
最近学了一下如何利用Fortran来对稀疏矩阵进行并行化计算，并且得到指定数目了一些本征态和本征值。
{:.info}
<!--more-->
# 矩阵产生
```fortran
! Author:YuXuanLi
! E-Mail:yxli406@gmail.com
! BHZ 模型，将Sz当作层自由度
    module pub
    implicit none
    integer xn,yn,N,len2
    parameter(xn = 30,yn = xn,N = xn*yn*8,len2 = xn*yn)
    complex,parameter::im = (0.0,1.0) 
    real,parameter::pi = 3.14159265359
    complex Ham(N,N) 
    integer bry2(4,len2)
    real m0,omega!,mu
    real tx,ty,del 
    real d0,dx,dy  
    real ax,ay,h0,tp 
!-------------------lapack parameter----------------------------------------
    integer::lda = N
    integer,parameter::lwmax=2*N+N**2
    real,allocatable::w(:)
    complex*8,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   
    integer lrwork    
    integer liwork   
    integer info
    end module pub
!========== PROGRAM START ==========================
    program sol
    use pub
    real t1,t2
    !================ Physics memory allocate =================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1 + 5*N + 2*N**2))
    allocate(iwork(3 + 5*N))
    !--------------------
    ! m0 = 1.
    ! tx = 0.4   
    ! ty = 1.3  
    ! ! mu = -0.3 
    ! ax = 0.4  
    ! ay = 1.3 
    m0 = 1.
    tx = 2.0   
    ty = 2.0  
    ! mu = -0.3 
    ax = 2.0
    ay = 2.0 
    h0 = 0. ! 层间耦合 
    tp = 0. !  inversion breaking
!=========== D(k) = D0 + D_x cos(k_x) + D_y cos(k_y) ============
    d0 = 0.5
    dx = 0.0
    dy = dx
    call cpu_time(t1)
    call matset(0.0)
    !call ldos(1)
    !call plot(1)
    call cpu_time(t2)
    open(20,file="time.dat")
    write(20,*)t2 - t1
    close(20)
    !call main1()
    stop
    end program sol
!============================================================
    subroutine main1()
    use pub
    integer i1
    i1 = 0
    do h0 = -1,1,0.1
        i1 = i1 + 1
        call matset(i1)
        call ldos(i1)
        call plot(i1)
    end do
    return
    end subroutine
!===========================Local Density of State=============================
    subroutine ldos(m2)
    use pub
    integer m,l1,k,l2,m2
    real s,omg
    character*20::str1,str2,str3,str4
    real,external::delta
    omg = 0.0
    str1 = "ldos-sc-"
    write(str2,"(I2.2)")m2
    str3 = ".dat"
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(12,file=str4) 
    do l1 = 1,yn
        do l2 = 1,xn
            k = l2 + (l1-1)*xn
            s = 0
            do m=1,N
                s = s + delta(w(m) - omg)*(abs(Ham(k , m))**2 + abs(Ham(k + len2 , m))**2 + abs(Ham(k + len2*2 , m))**2 + abs(Ham(k + len2*3 , m))**2)&
                + delta(w(m) - omg)*(abs(Ham(k + len2*4, m))**2 + abs(Ham(k + len2*5 , m))**2 + abs(Ham(k + len2*6 , m))**2 + abs(Ham(k + len2*7 , m))**2)
            end do
            write(12,*)l1,l2,s
        end do
        write(12,*)""
    end do
    close(12)
    return
    end subroutine ldos
!======================================================================
    real function delta(x)
    real x
    real::gamma = 0.001
    delta = 1.0/3.1415926535*gamma/(x*x + gamma*gamma)
    end function delta
!================= Topological insulator Hamiltonian ==============
    subroutine matset(input)
    use pub
    integer m,l,i,input 
    real mu
    Ham = 0
    call boundary()
    do m = 1,yn
        do l = 1,xn
            ! if(m<l)then
                i = (m - 1)*xn + l
                Ham(i , i) = m0 - mu 
                Ham(len2 + i , len2 + i) = -m0 - mu
                Ham(len2*2 + i , len2*2 + i) = m0 - mu  
                Ham(len2*3 + i , len2*3 + i) = -m0 - mu  
                Ham(len2*4 + i , len2*4 + i) = -m0 + mu  
                Ham(len2*5 + i , len2*5 + i) = m0 + mu  
                Ham(len2*6 + i , len2*6 + i) = -m0 + mu  
                Ham(len2*7 + i , len2*7 + i) = m0 + mu  
                !-----------------------------------------------
                Ham(i,len2*2 + i) = h0
                ham(len2 + i,len2*3 + i) = h0
                ham(len2*2 + i,i) = h0
                ham(len2*3 + i,len2 + i) = h0
                ham(len2*4 + i,len2*6 + i) = -h0
                ham(len2*5 + i,len2*7 + i) = -h0
                ham(len2*6 + i,len2*4 + i) = -h0
                ham(len2*7 + i,len2*5 + i) = -h0
                !-----------------------------------------------
                Ham(i,len2 + i) = tp
                ham(len2 + i, i) = tp
                ham(len2*2 + i,len2*3 + i) = tp
                ham(len2*3 + i,len2*2 + i) = tp
                ham(len2*4 + i,len2*5 + i) = -tp
                ham(len2*5 + i,len2*4 + i) = -tp
                ham(len2*6 + i,len2*7 + i) = -tp
                ham(len2*7 + i,len2*6 + i) = -tp
                !----------------------------------------------------------
                !(1,1)
                if (l.ne.xn)Ham(i,bry2(1,i)) = -tx/2
                if (l.ne.1)Ham(i,bry2(2,i)) = -tx/2
                if (m.ne.yn)Ham(i,bry2(3,i)) = -ty/2
                if (m.ne.1)Ham(i,bry2(4,i)) = -ty/2 
                !(2,2)
                if (l.ne.xn)Ham(len2 + i, len2 + bry2(1,i)) = tx/2
                if (l.ne.1)Ham(len2 + i, len2 + bry2(2,i)) = tx/2
                if (m.ne.yn)Ham(len2 + i, len2 + bry2(3,i)) = ty/2
                if (m.ne.1)Ham(len2 + i, len2 + bry2(4,i)) = ty/2
                !(3,3)
                if (l.ne.xn)Ham(len2*2 + i, len2*2 + bry2(1,i)) = -tx/2
                if (l.ne.1)Ham(len2*2 + i, len2*2 + bry2(2,i)) = -tx/2
                if (m.ne.yn)Ham(len2*2 + i, len2*2 + bry2(3,i)) = -ty/2
                if (m.ne.1)Ham(len2*2 + i, len2*2 + bry2(4,i)) = -ty/2
                !(4,4)
                if (l.ne.xn)Ham(len2*3 + i, len2*3 + bry2(1,i)) = tx/2
                if (l.ne.1)Ham(len2*3 + i, len2*3 + bry2(2,i)) = tx/2
                if (m.ne.yn)Ham(len2*3 + i, len2*3 + bry2(3,i)) = ty/2
                if (m.ne.1)Ham(len2*3 + i, len2*3 + bry2(4,i)) = ty/2
                !(5,5)
                if (l.ne.xn)Ham(len2*4 + i, len2*4 + bry2(1,i)) = tx/2
                if (l.ne.1)Ham(len2*4 + i, len2*4 + bry2(2,i)) = tx/2
                if (m.ne.yn)Ham(len2*4 + i, len2*4 + bry2(3,i)) = ty/2
                if (m.ne.1)Ham(len2*4 + i, len2*4 + bry2(4,i)) = ty/2
                !(6,6)
                if (l.ne.xn)Ham(len2*5 + i, len2*5 + bry2(1,i)) = -tx/2
                if (l.ne.1)Ham(len2*5 + i, len2*5 + bry2(2,i)) = -tx/2
                if (m.ne.yn)Ham(len2*5 + i, len2*5 + bry2(3,i)) = -ty/2
                if (m.ne.1)Ham(len2*5 + i, len2*5 + bry2(4,i)) = -ty/2
                !(7,7)
                if (l.ne.xn)Ham(len2*6 + i, len2*6 + bry2(1,i)) = tx/2
                if (l.ne.1)Ham(len2*6 + i, len2*6 + bry2(2,i)) = tx/2
                if (m.ne.yn)Ham(len2*6 + i, len2*6 + bry2(3,i)) = ty/2
                if (m.ne.1)Ham(len2*6 + i, len2*6 + bry2(4,i)) = ty/2
                !(8,8)
                if (l.ne.xn)Ham(len2*7 + i, len2*7 + bry2(1,i)) = -tx/2
                if (l.ne.1)Ham(len2*7 + i, len2*7 + bry2(2,i)) = -tx/2
                if (m.ne.yn)Ham(len2*7 + i, len2*7 + bry2(3,i)) = -ty/2
                if (m.ne.1)Ham(len2*7 + i, len2*7 + bry2(4,i)) = -ty/2
                !=========================================================
                !(1,2)
                if (l.ne.xn)Ham(i, len2 + bry2(1,i)) = -im*ax/2
                if (l.ne.1)Ham(i, len2 + bry2(2,i)) = im*ax/2 
                if (m.ne.yn)Ham(i, len2 + bry2(3,i)) = -ay/2
                if (m.ne.1)Ham(i, len2 + bry2(4,i)) = ay/2 
                !(2,1)
                if (l.ne.xn)Ham(len2 + i,bry2(1,i)) = -im*ax/2
                if (l.ne.1)Ham(len2 + i,bry2(2,i)) = im*ax/2 
                if (m.ne.yn)Ham(len2 + i,bry2(3,i)) = ay/2
                if (m.ne.1)Ham(len2 + i,bry2(4,i)) = -ay/2
                !(3,4)
                if (l.ne.xn)Ham(len2*2 + i, len2*3 + bry2(1,i)) = im*ax/2 
                if (l.ne.1)Ham(len2*2 + i, len2*3 + bry2(2,i)) = -im*ax/2 
                if (m.ne.yn)Ham(len2*2 + i, len2*3 + bry2(3,i)) = -ay/2
                if (m.ne.1)Ham(len2*2 + i, len2*3 + bry2(4,i)) = ay/2
                !(4,3)
                if (l.ne.xn)Ham(len2*3 + i, len2*2 + bry2(1,i)) = im*ax/2 
                if (l.ne.1)Ham(len2*3 + i, len2*2 + bry2(2,i)) = -im*ax/2 
                if (m.ne.yn)Ham(len2*3 + i, len2*2 + bry2(3,i)) = ay/2
                if (m.ne.1)Ham(len2*3 + i, len2*2 + bry2(4,i)) = -ay/2
                !(5,6)
                if (l.ne.xn)Ham(len2*4 + i, len2*5 + bry2(1,i)) = -im*ax/2 
                if (l.ne.1)Ham(len2*4 + i, len2*5 + bry2(2,i)) = im*ax/2 
                if (m.ne.yn)Ham(len2*4 + i, len2*5 + bry2(3,i)) = ay/2
                if (m.ne.1)Ham(len2*4 + i, len2*5 + bry2(4,i)) = -ay/2
                !(6,5)
                if (l.ne.xn)Ham(len2*5 + i, len2*4 + bry2(1,i)) = -im*ax/2 
                if (l.ne.1)Ham(len2*5 + i, len2*4 + bry2(2,i)) = im*ax/2 
                if (m.ne.yn)Ham(len2*5 + i, len2*4 + bry2(3,i)) = -ay/2
                if (m.ne.1)Ham(len2*5 + i, len2*4 + bry2(4,i)) = ay/2
                !(7.8)
                if (l.ne.xn)Ham(len2*6 + i, len2*7 + bry2(1,i)) = im*ax/2 
                if (l.ne.1)Ham(len2*6 + i, len2*7 + bry2(2,i)) = -im*ax/2 
                if (m.ne.yn)Ham(len2*6 + i, len2*7 + bry2(3,i)) = ay/2
                if (m.ne.1)Ham(len2*6 + i, len2*7 + bry2(4,i)) = -ay/2
                !(8,7)
                if (l.ne.xn)Ham(len2*7 + i, len2*6 + bry2(1,i)) = im*ax/2 
                if (l.ne.1)Ham(len2*7 + i, len2*6 + bry2(2,i)) = -im*ax/2 
                if (m.ne.yn)Ham(len2*7 + i, len2*6 + bry2(3,i)) = -ay/2
                if (m.ne.1)Ham(len2*7 + i, len2*6 + bry2(4,i)) = ay/2
                !============================================
                !(1,7)
                Ham(i, len2*6 + i) = -d0
                if (l.ne.xn)Ham(i, len2*6 + bry2(1,i)) = -dx/2
                if (l.ne.1)Ham(i, len2*6 + bry2(2,i)) = -dx/2
                if (m.ne.yn)Ham(i, len2*6 + bry2(3,i)) = -dy/2
                if (m.ne.1)Ham(i, len2*6 + bry2(4,i)) = -dy/2
                !(2,8)
                Ham(len2 + i, len2*7 + i) = -d0
                if (l.ne.xn)Ham(len2 + i, len2*7 + bry2(1,i)) = -dx/2
                if (l.ne.1)Ham(len2 + i, len2*7 + bry2(2,i)) = -dx/2
                if (m.ne.yn)Ham(len2 + i, len2*7 + bry2(3,i)) = -dy/2
                if (m.ne.1)Ham(len2 + i, len2*7 + bry2(4,i)) = -dy/2
                !(3,5)
                Ham(len2*2 + i, len2*4 + i) = d0
                if (l.ne.xn)Ham(len2*2 + i, len2*4 + bry2(1,i)) = dx/2
                if (l.ne.1)Ham(len2*2 + i, len2*4 + bry2(2,i)) = dx/2
                if (m.ne.yn)Ham(len2*2 + i, len2*4 + bry2(3,i)) = dy/2
                if (m.ne.1)Ham(len2*2 + i, len2*4 + bry2(4,i)) = dy/2
                !(4,6)
                Ham(len2*3 + i, len2*5 + i) = d0
                if (l.ne.xn)Ham(len2*3 + i, len2*5 + bry2(1,i)) = dx/2
                if (l.ne.1)Ham(len2*3 + i, len2*5 + bry2(2,i)) = dx/2
                if (m.ne.yn)Ham(len2*3 + i, len2*5 + bry2(3,i)) = dy/2
                if (m.ne.1)Ham(len2*3 + i, len2*5 + bry2(4,i)) = dy/2
                !(5,3)
                Ham(len2*4 + i,len2*2 + i) = d0
                if (l.ne.xn)Ham(len2*4 + i,len2*2 + bry2(1,i)) = dx/2
                if (l.ne.1)Ham(len2*4 + i,len2*2 + bry2(2,i)) = dx/2
                if (m.ne.yn)Ham(len2*4 + i,len2*2 + bry2(3,i)) = dy/2
                if (m.ne.1)Ham(len2*4 + i,len2*2 + bry2(4,i)) = dy/2
                !(6,4)
                Ham(len2*5 + i, len2*3 + i) = d0
                if (l.ne.xn)Ham(len2*5 + i, len2*3 + bry2(1,i)) = dx/2
                if (l.ne.1)Ham(len2*5 + i, len2*3 + bry2(2,i)) = dx/2
                if (m.ne.yn)Ham(len2*5 + i, len2*3 + bry2(3,i)) = dy/2
                if (m.ne.1)Ham(len2*5 + i, len2*3 + bry2(4,i)) = dy/2
                !(7,1)
                Ham(len2*6 + i, i) = -d0
                if (l.ne.xn)Ham(len2*6 + i, bry2(1,i)) = -dx/2
                if (l.ne.1)Ham(len2*6 + i, bry2(2,i)) = -dx/2
                if (m.ne.yn)Ham(len2*6 + i, bry2(3,i)) = -dy/2
                if (m.ne.1)Ham(len2*6 + i, bry2(4,i)) = -dy/2
                !(8,2)
                Ham(len2*7 + i, len2 + i) = -d0
                if (l.ne.xn)Ham(len2*7 + i, len2 + bry2(1,i)) = -dx/2
                if (l.ne.1)Ham(len2*7 + i, len2 + bry2(2,i)) = -dx/2
                if (m.ne.yn)Ham(len2*7 + i, len2 + bry2(3,i)) = -dy/2
                if (m.ne.1)Ham(len2*7 + i, len2 + bry2(4,i)) = -dy/2
            ! end if
        end do
    end do
    call isHermitian()
    call eigsol(input)
    open(10,file = "ham-re.dat")
    open(11,file = "ham-im.dat")
    write(10,999)real(Ham)
    write(11,999)aimag(Ham)
    close(10)
    close(11)
    ! do i1 = 1,N
    !     do i2 = 1,N
    !         write(10)
    !     end do
    ! end do
999 format(7200(7200f11.5/))
    return
    end subroutine matset
!===========================================================
    subroutine boundary()
    use pub
    integer i,ix,iy
    bry1 = 0
    do iy = 1,yn  
        do ix = 1,xn 
            i = (iy-1)*xn + ix
            bry2(1,i) = i + 1    
            if(ix.eq.xn)bry2(1,i) = bry2(1,i) - xn    
            bry2(2,i) = i - 1    
            if(ix.eq.1)bry2(2,i) = bry2(2,i) + xn     
            bry2(3,i) = i + xn   
            if(iy.eq.yn)bry2(3,i) = bry2(3,i) - len2  
            bry2(4,i)= i - xn    
            if(iy.eq.1)bry2(4,i) = bry2(4,i) + len2   
        enddo
    enddo
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
                write(*,*)"Hamiltonian is not Hermitian"
                stop
            end if
        end do
    end do
    close(16)
    return
    end subroutine isHermitian
!================= Hermitain Matrices solve ==============
    subroutine eigsol(input)
    use pub
    integer m,input
    character*20::str1,str2,str3,str4
    str1 = "eigval"
    str3 = ".dat"
    write(str2,"(I4.4)")input
    str4 = trim(str1)//trim(str2)//trim(str3)
    lwork = -1
    liwork = -1
    lrwork = -1
    call cheevd('V','U',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    lwork = min(2*N+N**2, int( work( 1 ) ) )
    lrwork = min(1+5*N+2*N**2, int( rwork( 1 ) ) )
    liwork = min(3+5*N, iwork( 1 ) )
    call cheevd('V','U',N,Ham,lda,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    if( info .GT. 0 ) then
        open(110,file="mes.dat",status="unknown")
        write(110,*)'The algorithm failed to compute eigenvalues.'
        close(110)
    end if
    open(120,file = str4) 
    do m = 1,N
        write(120,*)m,w(m)
    end do
    close(120)
    return
    end subroutine eigsol
!=======================================================================================
    subroutine plot(m3)
    use pub 
    character*20::str1,str2,str3,str4,str5,str6
    integer m3
    str1 = "ldos-sc-"
    write(str2,"(I2.2)")m3
    str3 = ".gnu"
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(20,file=str4)
    write(20,*)'set encoding iso_8859_1'
    write(20,*)'#set terminal  postscript enhanced color'
    write(20,*)"#set output 'arc_r.eps'"
    write(20,*)'#set terminal pngcairo truecolor enhanced  font ",50" size 1920, 1680'
    write(20,*)'set terminal png truecolor enhanced font ",50" size 1920, 1680'
    write(20,*)"set output 'ldos-sc-"//trim(str2)//".png'"
    write(20,*)'#set palette defined ( -10 "#194eff", 0 "white", 10 "red" )'
    write(20,*)'set palette defined ( -10 "blue", 0 "white", 10 "red" )'
    write(20,*)'#set palette rgbformulae 33,13,10'
    write(20,*)'unset ztics'
    write(20,*)'unset key'
    write(20,*)'set pm3d'
    write(20,*)'set border lw 6'
    write(20,*)'set size ratio -1'
    write(20,*)'set view map'
    write(20,*)'set xtics'
    write(20,*)'set ytics'
    write(str5,'(f5.2)')h0 ! 将一个小数转成字符串
    write(20,*)"set title 'h0 ="//trim(str5)//"'"
    write(20,*)'#set xlabel "K_1 (1/{\305})"'
    write(20,*)'set xlabel "X"'
    write(20,*)'#set ylabel "K_2 (1/{\305})"'
    write(20,*)'set ylabel "Y"'
    write(20,*)'set ylabel offset 1, 0'
    write(20,*)'set colorbox'
    write(20,*)'set xrange [1: 30]'
    write(20,*)'set yrange [1: 30]'
    write(20,*)'set pm3d interpolate 4,4'
    write(20,*)"splot 'ldos-sc-"//trim(str2)//".dat' u 1:2:3 w pm3d"
    close(20)
    !-------------------
    str1 = "val-sc-"
    write(str2,"(I4.4)")m3
    str3 = ".gnu"
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(21,file=str4)
    write(21,*)'set encoding iso_8859_1'
    write(21,*)'#set terminal  postscript enhanced color'
    write(21,*)"#set output 'arc_r.eps'"
    write(21,*)'#set terminal pngcairo truecolor enhanced  font ",50" size 1920, 1680'
    write(21,*)'set terminal png truecolor enhanced font ",50" size 1920, 1680'
    write(21,*)"set output 'val-sc-"//trim(str2)//".png'"
    write(21,*)"unset key"
    write(str5,'(f5.2)')h0 ! 将一个小数转成字符串
    write(21,*)"set title 'h0 ="//trim(str5)//"'"
    write(21,*)'set xrange [3540:3660]'
    write(21,*)'set yrange [-0.5:0.5]'
    write(21,*)"plot 'eigval"//trim(str2)//".dat' u 1:2 ps 3.0 pt 7"
    close(21)
    return 
    end subroutine plot
```
其中最主要的就是构建出一个哈密顿量，然后将其实部和虚部分别写出来
```fortran
  open(10,file = "ham-re.dat")
  open(11,file = "ham-im.dat")
  write(10,999)real(Ham)
  write(11,999)aimag(Ham)
  close(10)
  close(11)
999 format(7200(7200f11.5/))
```
# 稀疏矩阵对角化
```fortran
      PROGRAM SAMPLE_PZHEEVX_call
      INTEGER            LWORK, MAXN, LIWORK
      doUBLE PRECISION   ZERO
      PARAMETER          ( LWORK = 100000000, MAXN = 46000,
     $                   LIWORK = 500000,
     $                   ZERO = 0.0E+0 )
      INTEGER            LDA
      doUBLE PRECISION   MONE
      INTEGER            MAXPROCS
      PARAMETER          ( LDA = MAXN, MONE = -1.0D+0, MAXPROCS = 512 )
      INTEGER            LRWORK
      PARAMETER          ( LRWORK = 100000000 )
*     ..
*     .. Local Scalars ..
      INTEGER            CONTEXT, I, IAM, INFO, M, MYCOL, MYROW, N, NB,
     $                   NPCOL, NPROCS, NPROW, NZ, vallen
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( 50 ), DESCZ( 50 ),
     $                   ICLUSTR( MAXPROCS*2 ), ifAIL( MAXN ),
     $                   IWORK( LIWORK )
      doUBLE PRECISION   GAP( MAXPROCS ), RWORK( LRWORK ), W( MAXN )
      COMPLEX*16         A( LDA, LDA ), WORK( LWORK ), Z( LDA, LDA )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP, DESCINIT, PZHEEVX, PZLAMODHILB,
     $                   PZLAPRNT
*     ..
*     .. Executable Statements ..
*
*
*     Set up the problem
*
      N = 7200 ! 哈密顿量矩阵的维度
      NB = 32
      NPROW = 4
      NPCOL = 4
      vallen = 40 ! 求解本征态的数量
*
*
*     Initialize the BLACS
*
      call BLACS_PINFO( IAM, NPROCS )
      if( ( NPROCS.LT.1 ) ) then
         call BLACS_SETUP( IAM, NPROW*NPCOL )
      end if
*
*
*     Initialize a single BLACS context
*
      call BLACS_GET( -1, 0, CONTEXT )
      call BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      call BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Bail out if this process is not a part of this context.
*
      if( MYROW.EQ.-1 )
     $   GO TO 20
*
*
*     These are basic array descriptors
*
      call DESCINIT( DESCA, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
      call DESCINIT( DESCZ, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
*
*     Build a matrix that you can create with
*     a one line matlab command:  hilb(n) + diag([1:-1/n:1/n])
*
      call PZLAMODHILB( N, A, 1, 1, DESCA, INFO )
*
*
*     Uncomment this line to see the matrix printed out.
*
*      call PZLAPRNT( N, N, A, 1, 1, DESCA, 0, 0, 'A', 6, WORK )
*
*
*     Ask PZHEEVX to compute the entire eigendecomposition
*
      call PZHEEVX( 'V', 'I', 'U', N, A, 1, 1, DESCA, ZERO, ZERO, 
     $              N/2 - vallen, N/2 + vallen, 
     $              MONE, M, NZ, W, MONE, Z, 1, 1, DESCZ, WORK,
     $              LWORK, RWORK, LRWORK, IWORK, LIWORK, ifAIL, ICLUSTR,
     $              GAP, INFO )
*
*
*     Print out the eigenvectors
      call PZLAPRNT( N, 21, Z, 1, 1, DESCZ, 0, 0, 'Z', 6, WORK )
      call BLACS_GRIDEXIT( CONTEXT )
*
   20 CONTINUE
      call BLACS_EXIT( 0 )
*     Uncomment this line on SUN systems to avoid the useless print out
*
*      call IEEE_FLAGS( 'clear', 'exception', 'underflow', '')
 9999 FORMAT( 'W=diag([', 4f16.12, ']);' )
      STOP
      end
!==============================================================================
      SUBROUTINE PZLAMODHILB( N, A, IA, JA, DESCA, INFO )
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      doUBLE PRECISION   ONE, HR(N), HI(N)
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      COMPLEX*16         A( * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PZELSET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX
      if( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DT_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )return
*
      INFO = 0
*
      call BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
*
      if( IA.NE.1 ) then
         INFO = -3
      else if( JA.NE.1 ) then
         INFO = -4
      end if
*
      Open(15, File="ham-re.dat", status='old')
      Open(12, File="ham-im.dat", status='old')
      do 20 J = 1, N
         read(15, *) (HR(JJ), JJ = 1,N)
         read(12, *) (HI(JJ), JJ = 1,N)
         do 10 I = 1, N
               call PZELSET( A, I, J, DESCA, dcmplx(HR(I), HI(I)))
   10    CONTINUE
   20 CONTINUE
      CLOSE(15)
      CLOSE(12)
      return
      end 
```
这里使用的是ScalPack的函数，只需要改动一下构建哈密顿量的这部分就可以了
```fortran
 Open(15, File="ham-re.dat", status='old')
      Open(12, File="ham-im.dat", status='old')
      do 20 J = 1, N
         read(15, *) (HR(JJ), JJ = 1,N)
         read(12, *) (HI(JJ), JJ = 1,N)
         do 10 I = 1, N
               call PZELSET( A, I, J, DESCA, dcmplx(HR(I), HI(I)))
   10    CONTINUE
   20 CONTINUE
```
在对角化的过程中，主要函数为
```fortran
call PZHEEVX( 'V', 'I', 'U', N, A, 1, 1, DESCA, ZERO, ZERO, 
     $              N/2 - vallen, N/2 + vallen, 
     $              MONE, M, NZ, W, MONE, Z, 1, 1, DESCZ, WORK,
     $              LWORK, RWORK, LRWORK, IWORK, LIWORK, ifAIL, ICLUSTR,
     $              GAP, INFO )
*
```
这里就是要的
> N/2 - vallen, N/2 + vallen

个本征值以及对应的本征态，这里计算得到的本征态会自动打印出来，还需要再处理一下就可以的到想要的本征态的
```shell
cat lyx.o356601 |grep 'Z (' >out
awk '{print $5,$9}' out>out2
```
这里假设程序计算的输出都存储到了`lyx.o356601`这个文件中，那么对这个文件进行上面的操作之后，就可以在`out2`中的到对应本征值的本征矢量。