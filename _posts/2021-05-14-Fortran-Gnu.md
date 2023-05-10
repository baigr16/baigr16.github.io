---
title: Fortran + Gnu 批量计算
tags: Fortran Plot Code Shell
layout: article
license: true
toc: true
key: a20210514
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
这篇博客主要介绍如何利用Fortran与Gunplot结合进行批量计算,主要针对模型研究时连续改变某个参数,批量的绘制想要结果图,这可以节省很多的时间,主要是在利用Fortran计算的时候,将相应参数绘图的**gnuplot**绘图代码写入相应文件,之后进行快速执行.
{:.info}
<!--more-->
做模型研究的时候,通常会调节一个模型的很多参数来研究这些参数变化时某些量是如何演化的,这就涉及到批量计算绘图的问题,我是习惯用Fortran了,因为计算速度比较快,而且自己也比较熟悉了,结合**gnuplot**可以快速的将不同数据对应的结果进行作图.这里的有些内容可以参考[做数值计算好用的软件及杂项整理](https://yxli8023.github.io/2020/09/16/introduction.html)这篇中的内容,废话不多说,直接上代码进行解释.

# Fortran + Gnuplot 批量数据输出绘图
```fortran
!   Anticommute mass term along y open gap
!   tau*spin*orbit
!   Rotation 45 degree calculate edge states
!   kx--->(Kx + Ky)    ky----->(Kx - Ky)
    module pub
    implicit none
    integer yn,kn,ne,N,enn,hn
    real eng   ! energy
    parameter(yn = 50,hn =  8,kn = 50, ne = 50,N = yn*hn)
    real,parameter::pi = 3.1415926535
    complex,parameter::im = (0.,1.0)  
    complex Ham(N,N) 
    !--------------------------------
    real m0,mu,lamx,lamy,gamma,tx,ty,dx,dy,d0
    complex g1(hn,hn) ,g4(hn,hn),g2(hn,hn),g3(hn,hn),g5(hn,hn)
    complex am1(hn,hn),am2(hn,hn),am3(hn,hn),am4(hn,hn)
    complex am5(hn,hn),am6(hn,hn),am7(hn,hn),am8(hn,hn),am(hn,hn)
    real mm0,mmx,mmy
    !---------------------------------
    integer::lda = N
    integer,parameter::lwmax = 2*N + N**2
    real,allocatable::w(:)
    complex,allocatable::work(:)
    real,allocatable::rwork(:)
    integer,allocatable::iwork(:)
    integer lwork   
    integer lrwork    
    integer liwork   
    integer info
    end module pub
!================= PROGRAM START ============================
    program sol
    use pub
!====空间申请==================
    allocate(w(N))
    allocate(work(lwmax))
    allocate(rwork(1+5*N+2*N**2))
    allocate(iwork(3+5*N))
!======parameter value setting =====
    m0 = 1.5
    ! mu = 0.3
    tx = 1.0
    ty = 1.0
    lamx = 1.0
    lamy = 1.0
    ! dx = 0.5
    ! dy = -dx
    mm0 = 0.2
    ! mmx = 0.5
    ! mmy = -mmx
    call Pauli()
    call main1()
    ! am = am1
    ! call cylinder(1)
    ! call plot(1)
    stop 
    end program sol
!====================================================================================
    subroutine main1()
    use pub
    am = am1
    call cylinder(1)
    call plot(1)

    am = am2
    call cylinder(2)
    call plot(2)

    am = am3
    call cylinder(3)
    call plot(3)

    am = am4
    call cylinder(4)
    call plot(4)

    am = am5
    call cylinder(5)
    call plot(5)

    am = am6
    call cylinder(6)
    call plot(6)

    am = am7
    call cylinder(7)
    call plot(7)

    am = am8
    call cylinder(8)
    call plot(8)

    end subroutine main1
!=============================================================================================
    subroutine cylinder(m3)
    !  Calculate edge spectrum function
    use pub
    integer m1,m2,m3
    real kx,ky
    character*20::str1,str2,str3,str4,str5,str6
    str1 = "oy45-sc-ym"
    str2 = "ox45-sc-ym"
    write(str3,"(I2.2)")m3
    str4 = ".dat"
    str5 = trim(str1)//trim(str3)//trim(str4)
    str6 = trim(str2)//trim(str3)//trim(str4)
    open(30,file=str5)
    !-------------------------------------------------
    !   y-direction is open
    do m1 = -kn,kn
        kx = pi*m1/kn
        call openy(kx)
        write(30,999)kx/pi,(w(i),i = 1,N)
    end do
    close(30)
    !--------------------------------------------------
    !  x-direction is open
    open(31,file=str6)
    do m1 = -kn,kn
        ky = pi*m1/kn
        call openx(ky)
        write(31,999)ky/pi,(w(i),i=1,N)
    end do
    close(31)
999 format(802f11.5)
    return
    end subroutine cylinder
!======================== Pauli Matrix driect product============================
    subroutine Pauli()
    use pub
    !---------Kinetic energy
    g1(1,1) = 1
    g1(2,2) = -1
    g1(3,3) = 1
    g1(4,4) = -1
    g1(5,5) = -1
    g1(6,6) = 1
    g1(7,7) = -1
    g1(8,8) = 1
    !----------SOC-x
    g2(1,2) = 1
    g2(2,1) = 1
    g2(3,4) = -1
    g2(4,3) = -1
    g2(5,6) = 1
    g2(6,5) = 1
    g2(7,8) = -1
    g2(8,7) = -1
    !---------SOC-y
    g3(1,2) = -im
    g3(2,1) = im
    g3(3,4) = -im
    g3(4,3) = im
    g3(5,6) = im 
    g3(6,5) = -im
    g3(7,8) = im
    g3(8,7) = -im
    !-------------------SC pairing
    g4(1,7) = -1
    g4(2,8) = -1
    g4(3,5) = 1
    g4(4,6) = 1
    g4(5,3) = 1
    g4(6,4) = 1
    g4(7,1) = -1
    g4(8,2) = -1
    !------------------- Chemical potential
    g5(1,1) = 1
    g5(2,2) = 1
    g5(3,3) = 1
    g5(4,4) = 1
    g5(5,5) = -1
    g5(6,6) = -1
    g5(7,7) = -1
    g5(8,8) = -1
    !-------------------Anti Mass-------------------------------
    !i,y,i
    am1(1,3) = -im
    am1(2,4) = -im
    am1(3,1) = im
    am1(4,2) = im
    am1(5,7) = -im
    am1(6,8) = -im
    am1(7,5) = im
    am1(8,6) = im
    !y,x,x
    am2(1,8) = -im
    am2(2,7) = -im
    am2(3,6) = -im
    am2(4,5) = -im
    am2(5,4) = im
    am2(6,3) = im
    am2(7,2) = im
    am2(8,1) = im
    !i,x,i
    am3(1,3) = 1
    am3(2,4) = 1
    am3(3,1) = 1
    am3(4,2) = 1
    am3(5,7) = 1
    am3(6,8) = 1
    am3(7,5) = 1
    am3(8,6) = 1
    !z,y,i
    am4(1,3) = -im
    am4(2,4) = -im
    am4(3,1) = im
    am4(4,2) = im
    am4(5,7) = im
    am4(6,8) = im
    am4(7,5) = -im
    am4(8,6) = -im
    !x,y,x
    am5(1,8) = -im
    am5(2,7) = -im
    am5(3,6) = im
    am5(4,5) = im
    am5(5,4) = -im
    am5(6,3) = -im
    am5(7,2) = im
    am5(8,1) = im
    !z,x,i
    am6(1,3) = 1
    am6(2,4) = 1
    am6(3,1) = 1
    am6(4,2) = 1
    am6(5,7) = -1
    am6(6,8) = -1
    am6(7,5) = -1
    am6(8,6) = -1
    !x,x,x
    am7(1,8) = 1
    am7(2,7) = 1
    am7(3,6) = 1
    am7(4,5) = 1
    am7(5,4) = 1
    am7(6,3) = 1
    am7(7,2) = 1
    am7(8,1) = 1
    !y,y,x
    am8(1,8) = -1
    am8(2,7) = -1
    am8(3,6) = 1
    am8(4,5) = 1
    am8(5,4) = 1
    am8(6,3) = 1
    am8(7,2) = -1
    am8(8,1) = -1
    return
    end subroutine Pauli
!==========================================================
    subroutine openx(ky)
    use pub
    real ky
    integer m,l,k
    Ham = 0
!========== Positive energy ========
    do k = 0,yn-1
        if (k == 0) then ! Only right block in first line
            do m = 1,hn
                do l = 1,hn
                    Ham(m,l) = m0*g1(m,l)+ mm0*am(m,l)

                    Ham(m,l + hn) = (tx/2.0*cos(ky) - tx/(2*im)*sin(ky) + ty/2.0*cos(ky) + ty/(2*im)*sin(ky))*g1(m,l) +&
                                    (lamx/(2*im)*cos(ky) + lamx/2.0*sin(ky))*g2(m,l) + (lamy/(2*im)*cos(ky) - lamy/2.0*sin(ky))*g3(m,l)&
                                    + (mmx/2.0*cos(ky) + mmx/(2*im)*sin(ky) - mmy/2.0*cos(ky) - mmy/(2*im)*sin(ky))*am(m,l)
                end do
            end do
        elseif ( k==yn-1 ) then ! Only left block in last line
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) = m0*g1(m,l)+ mm0*am(m,l)

                    Ham(k*hn + m,k*hn + l - hn) = (tx/2.0*cos(ky) + tx/(2*im)*sin(ky) + ty/2.0*cos(ky) - ty/(2*im)*sin(ky))*g1(m,l) +&
                    (-lamx/(2*im)*cos(ky) + lamx/2.0*sin(ky))*g2(m,l) + (-lamy/(2*im)*cos(ky) - lamy/2.0*sin(ky))*g3(m,l)&
                    + (mmx/2.0*cos(ky) - mmx/(2*im)*sin(ky) - mmy/2.0*cos(ky) + mmy/(2*im)*sin(ky))*am(m,l)
                end do
            end do
        else
            do m = 1,hn
                do l = 1,hn ! k start from 1,matrix block from 2th row
                    Ham(k*hn + m,k*hn + l) = m0*g1(m,l)+ mm0*am(m,l)

                    Ham(k*hn + m,k*hn + l + hn) = (tx/2.0*cos(ky) - tx/(2*im)*sin(ky) + ty/2.0*cos(ky) + ty/(2*im)*sin(ky))*g1(m,l) +&
                    (lamx/(2*im)*cos(ky) + lamx/2.0*sin(ky))*g2(m,l) + (lamy/(2*im)*cos(ky) - lamy/2.0*sin(ky))*g3(m,l)&
                    + (mmx/2.0*cos(ky) + mmx/(2*im)*sin(ky) - mmy/2.0*cos(ky) - mmy/(2*im)*sin(ky))*am(m,l)
                    Ham(k*hn + m,k*hn + l - hn)  = (tx/2.0*cos(ky) + tx/(2*im)*sin(ky) + ty/2.0*cos(ky) - ty/(2*im)*sin(ky))*g1(m,l) +&
                    (-lamx/(2*im)*cos(ky) + lamx/2.0*sin(ky))*g2(m,l) + (-lamy/(2*im)*cos(ky) - lamy/2.0*sin(ky))*g3(m,l)&
                    + (mmx/2.0*cos(ky) - mmx/(2*im)*sin(ky) - mmy/2.0*cos(ky) + mmy/(2*im)*sin(ky))*am(m,l)
                end do
            end do
        end if
    end do
    !------------------------
    call isHermitian()
    call eigsol()
    return
    end subroutine openx
!==========================================================
    subroutine openy(kx)
    use pub
    real kx
    integer m,l,k
    Ham = 0
!========== Positive energy ========
    do k = 0,yn-1
        if (k == 0) then ! Only right block in first line
            do m = 1,hn
                do l = 1,hn
                    Ham(m,l) =  m0*g1(m,l) + mm0*am(m,l)

                    Ham(m,l + hn) = (tx/2.0*cos(kx) - tx/(2*im)*sin(kx) + ty/2.0*cos(kx) + ty/(2*im)*sin(kx))*g1(m,l)+&
                                    (lamx/2.0*sin(kx) + lamx/(2*im)*cos(kx))*g2(m,l) + (lamy/2.0*sin(kx) - lamy/(2*im)*cos(kx))*g3(m,l)&
                                    + (mmx/2.0*cos(kx) + mmx/(2*im)*sin(kx) - mmy/2.0*cos(kx) - mmy/(2*im)*sin(kx))*am(m,l)
                end do
            end do
        elseif ( k==yn-1 ) then ! Only left block in last line
            do m = 1,hn
                do l = 1,hn
                    Ham(k*hn + m,k*hn + l) =  m0*g1(m,l)+ mm0*am(m,l)

                    Ham(k*hn + m,k*hn + l - hn) = (tx/2.0*cos(kx) + tx/(2*im)*sin(kx) + ty/2.0*cos(kx) - ty/(2*im)*sin(kx))*g1(m,l)+&
                    (lamx/2.0*sin(kx) - lamx/(2*im)*cos(kx))*g2(m,l) + (lamy/2.0*sin(kx) + lamy/(2*im)*cos(kx))*g3(m,l)&
                    + (mmx/2.0*cos(kx) - mmx/(2*im)*sin(kx) - mmy/2.0*cos(kx) + mmy/(2*im)*sin(kx))*am(m,l)
                end do
            end do
        else
            do m = 1,hn
                do l = 1,hn ! k start from 1,matrix block from 2th row
                    Ham(k*hn + m,k*hn + l) =  m0*g1(m,l)+ mm0*am(m,l)

                    Ham(k*hn + m,k*hn + l + hn) = (tx/2.0*cos(kx) - tx/(2*im)*sin(kx) + ty/2.0*cos(kx) + ty/(2*im)*sin(kx))*g1(m,l)+&
                    (lamx/2.0*sin(kx) + lamx/(2*im)*cos(kx))*g2(m,l) + (lamy/2.0*sin(kx) - lamy/(2*im)*cos(kx))*g3(m,l)&
                    + (mmx/2.0*cos(kx) + mmx/(2*im)*sin(kx) - mmy/2.0*cos(kx) - mmy/(2*im)*sin(kx))*am(m,l)
                    Ham(k*hn + m,k*hn + l - hn) = (tx/2.0*cos(kx) + tx/(2*im)*sin(kx) + ty/2.0*cos(kx) - ty/(2*im)*sin(kx))*g1(m,l)+&
                    (lamx/2.0*sin(kx) - lamx/(2*im)*cos(kx))*g2(m,l) + (lamy/2.0*sin(kx) + lamy/(2*im)*cos(kx))*g3(m,l)&
                    + (mmx/2.0*cos(kx) - mmx/(2*im)*sin(kx) - mmy/2.0*cos(kx) + mmy/(2*im)*sin(kx))*am(m,l)
                end do
            end do
        end if
    end do
    !---------------------------------
    call isHermitian()
    call eigsol()
    return
    end subroutine openy
!============================================================
    subroutine isHermitian()
    use pub
    integer i,j
    do i = 1,N
        do j = 1,N
            if (Ham(i,j) .ne. conjg(Ham(j,i)))then
                open(160,file = 'hermitian.dat')
                write(160,*)i,j
                write(160,*)Ham(i,j)
                write(160,*)Ham(j,i)
                write(160,*)"===================="
                write(*,*)"Hamiltonian is not hermitian"
                stop
            end if
        end do
    end do
    close(160)
    return
    end subroutine isHermitian
!================= Hermitain Matrices solve ==============
    subroutine eigsol(input)
    use pub
    integer m
    character*20::str1,str2,str3,str4
    str1 = "eigval_xmass"
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
    return
    end subroutine eigsol
!==========================================================================
    subroutine plot(m3)
    use pub
    character*20::str1,str2,str3,str4,str5,str6
    integer m3
    str1 = "diag-ym"
    write(str2,"(I2.2)")m3
    str3 = ".gnu"
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(20,file=str4)
    write(20,*)'set terminal png truecolor enhanced font ",50" size 3000, 1500'
    write(20,*)"set output 'diag"//trim(str2)//"-cy-ym.png'"
    write(20,*)"set size 1, 1"
    write(20,*)"set multiplot layout 1, 2"
    write(20,*)"unset key"
    write(20,*)"set ytics 1.5 "
    write(20,*)"set xtics 0.5"
    write(20,*)"set xtics offset 0, 0.0"
    write(20,*)'set xtics format "%4.1f" nomirror out '
    write(20,*)'set ytics format "%4.1f"'
    write(20,*)'set ylabel "E"'
    write(20,*)"set ylabel offset 0.5, 0.0 "
    write(20,*)"#set xlabel offset 0, -1.0 "
    write(20,*)"set xrange [-1:1]"
    write(20,*)"set yrange [-3:3]"
    write(20,*)"set xlabel 'k_y' "
    write(20,*)"plot for [i=2:400] 'ox45-sc-ym"//trim(str2)//".dat' u 1:i w l lw 5 lc 'blue'"
    write(20,*)"set xlabel 'k_x' "
    write(20,*)"plot for [i=2:400] 'oy45-sc-ym"//trim(str2)//".dat' u 1:i w l lw 5 lc 'blue'"
    write(20,*)"unset multiplot"
    close(20)
    end subroutine plot
```

上面代码主要的内容是
```fortran
  subroutine cylinder(m3)
    !  Calculate edge spectrum function
    use pub
    integer m1,m2,m3
    real kx,ky
    character*20::str1,str2,str3,str4,str5,str6
    str1 = "oy45-sc-ym"
    str2 = "ox45-sc-ym"
    write(str3,"(I2.2)")m3
    str4 = ".dat"
    str5 = trim(str1)//trim(str3)//trim(str4)
    str6 = trim(str2)//trim(str3)//trim(str4)
    open(30,file=str5)
    !-------------------------------------------------
    !   y-direction is open
    do m1 = -kn,kn
        kx = pi*m1/kn
        call openy(kx)
        write(30,999)kx/pi,(w(i),i = 1,N)
    end do
    close(30)
    !--------------------------------------------------
    !  x-direction is open
    open(31,file=str6)
    do m1 = -kn,kn
        ky = pi*m1/kn
        call openx(ky)
        write(31,999)ky/pi,(w(i),i=1,N)
    end do
    close(31)
999 format(802f11.5)
    return
    end subroutine cylinder
```
通过下面的代码可以组合出随着**m3**变化的字符串,将其作为文件名,这样就会对不同的输出,产生相对应的数据文件
```fortran
    str1 = "oy45-sc-ym"
    str2 = "ox45-sc-ym"
    write(str3,"(I2.2)")m3
    str4 = ".dat"
    str5 = trim(str1)//trim(str3)//trim(str4)
    str6 = trim(str2)//trim(str3)//trim(str4)
```
而在绘图部分
```fortran
    subroutine plot(m3)
    use pub
    character*20::str1,str2,str3,str4,str5,str6
    integer m3
    str1 = "diag-ym"
    write(str2,"(I2.2)")m3
    str3 = ".gnu"
    str4 = trim(str1)//trim(str2)//trim(str3)
    open(20,file=str4)
    write(20,*)'set terminal png truecolor enhanced font ",50" size 3000, 1500'
    write(20,*)"set output 'diag"//trim(str2)//"-cy-ym.png'"
    write(20,*)"set size 1, 1"
    write(20,*)"set multiplot layout 1, 2"
    write(20,*)"unset key"
    write(20,*)"set ytics 1.5 "
    write(20,*)"set xtics 0.5"
    write(20,*)"set xtics offset 0, 0.0"
    write(20,*)'set xtics format "%4.1f" nomirror out '
    write(20,*)'set ytics format "%4.1f"'
    write(20,*)'set ylabel "E"'
    write(20,*)"set ylabel offset 0.5, 0.0 "
    write(20,*)"#set xlabel offset 0, -1.0 "
    write(20,*)"set xrange [-1:1]"
    write(20,*)"set yrange [-3:3]"
    write(20,*)"set xlabel 'k_y' "
    write(20,*)"plot for [i=2:400] 'ox45-sc-ym"//trim(str2)//".dat' u 1:i w l lw 5 lc 'blue'"
    write(20,*)"set xlabel 'k_x' "
    write(20,*)"plot for [i=2:400] 'oy45-sc-ym"//trim(str2)//".dat' u 1:i w l lw 5 lc 'blue'"
    write(20,*)"unset multiplot"
    close(20)
    end subroutine plot
```
同样利用相同的方法来产生相对应数据的**gnuplot**作图文件.

这里再提供一个编译执行fortran文件的shell脚本
```shell
#!/bin/sh  
#============ get the file name ===========  
rm *.gnu *.png no* *.dat *.out 1>/dev/null 2>/dev/null &
Folder_A=$(pwd) 
for file_a in ${Folder_A}/*.f90
do 
	temp_file=`basename $file_a  .f90` 
	ifort -mkl -O3 -CB $file_a -o $temp_file.out 
    nohup ./$temp_file.out 1>/dev/null 2>/dev/null &
done
```
将这个脚本文件放在放置在和需要编译的fortran文件相同的文件夹中,这个脚本会自动搜寻当前文件夹中所有后缀名为**.f90**的文件,然后编译并执行,这里采用的**ifort**进行编译的,而且因为对角化用到了mkl库中的**cheevd**这个函数,所有编译选项有-mkl这个参数.
{:.success}

将上面的脚本命名为**run1.sh**然后执行
```shell
sh run1.sh
```
结果如下图所示

![png](/assets/images/Fortran/for1.png)

可以看到这里得到了数据结果**filename.dat**与其相应的一系列绘图脚本**filename.gnu**,下面批量执行所有的**filename.gnu**绘图脚本
```shell
#!/bin/sh  
#============ get the file name ===========  
Folder_A=$(pwd) 
for file_a in ${Folder_A}/*.gnu
do 
	gnuplot $file_a  
done
```
同样的,这个脚本的作用是寻找该文件夹下面所有后缀为**.gnu**的文件,执行绘图程序**gnuplot filename.gnu**,结果如下

![png](/assets/images/Fortran/for2.png)

可以看到批量执行之后会有一系列的结果图.

除了上面对单个文件的操作,有时候如果是对单个文件夹中的多个文件操作,上面的操作同样是可性的,因为编译执行的脚本是对后缀名进行识别的,所以如果有多个fortran文件,那么它们都会被编译执行,而绘图的脚本也是同样的道理.
{:.warning}

# 文件夹递归操作
除了对单个文件夹中的多个文件进行编译执行操作,有时候在计算时,可能会有很多个文件夹,而每个文件夹中又有多个文件需要进行编译执行,这时候上面的脚本就不能发挥作用了,需要对文件夹进行递归操作,也就是说要找到当前文件夹下面包含的每个文件夹,进行该文件夹然后再执行上面的脚本,比如当前文件夹下面还有一些文件夹

![png](/assets/images/Fortran/for3.png)


```shell
#!/bin/bash 
# 递归搜寻文件夹下面所有的.f90或者.f后缀结尾的文件,并利用ifort编译该文件,然后执行
function getdir(){
    for element in `ls $1`
      do
        dir_or_file=$1"/"$element
    if [ -d $dir_or_file ]
      then
        getdir $dir_or_file
      else  # 下面的全是文件
	  	if [ "${dir_or_file##*.}"x = "f90"x ]||[ "${dir_or_file##*.}"x = "f"x ];then	# 筛选处特定后缀的文件
    		dir_name=`dirname $dir_or_file` # 读取目录
			file_name=`basename $dir_or_file .f90` # 读取以f90结尾的文件名
			out_file_name="$dir_name/$file_name"  # 定义编号成功的文件名称
			ifort -mkl $dir_or_file -o $out_file_name.out  # 开始编译fortran文件,编译后文件名以out结尾
			dir1=`dirname $out_file_name`
			#echo $dir1
			cd $dir1  # 切换到具体的文件夹
			sh run1.sh & # 执行当前文件夹下面的run1.sh脚本
			#./$file_name.out 1>mes 2>bad &  # 执行该文件夹下面编译好的文件
			# ./$out_file_name.out 1>mes 2>bad &
			rm $out_file_name.out
		fi
        #temp_file=`basename $dir_or_file  .f90` #将文件名后缀删除
        #ifort -mkl $dir_or_file -o $temp_file.out  # 编译后文件名以out结尾
        #echo $dir_or_file     # 这里的变量时完整的路径名
    fi
done
}
 
#root_dir="/home/yxli/te"
root_dir=`pwd`
getdir $root_dir
```
代码的主要部分是
```shell
     $dir1  # 切换到具体的文件夹
			sh run1.sh & # 执行当前文件夹下面的run1.sh脚本
			#./$file_name.out 1>mes 2>bad &  # 执行该文件夹下面编译好的文件
			# ./$out_file_name.out 1>mes 2>bad &
			rm $out_file_name.out
```
我在这里使用了**sh run1.sh**这个命令,这样的好处是会对cd切换到的每个文件夹中均编译执行所有**.f90**的文件(这是脚本**run1.sh**的功能),同样也可以使用
```shell
      #./$file_name.out 1>mes 2>bad &  # 执行该文件夹下面编译好的文件
			# ./$out_file_name.out 1>mes 2>bad &
```
这个命令,我这里是将它进行了注释,它也是会编译并执行所有后缀为**.f90**的文件,而关于后缀名的选择,是在
```shell
if [ "${dir_or_file##*.}"x = "f90"x ]||[ "${dir_or_file##*.}"x = "f"x ];then	# 筛选处特定后缀的文件
```
这里进行筛选的,所以如果是想要编译执行其他类型程序的文件,可以对这个后缀名进行修改并调整对应的编译器即可.

这个脚本的主要思路就是,遍历当前文件夹下的所有文件夹,而每个文件夹中又有前面提及到的**run1.sh**这个脚本,在循环遍历所有的文件夹过程中,执行每个文件夹下面的**run1.sh**脚本,这样就相当于对所有的文件夹都启动的编译执行程序,由于这是一个递归的过程,就算当前文件夹下面包含的文件夹深度不止一层,这个脚本同样可以使用.
{:.success}

在将fortran程序递归的编译执行完成之后,每个文件夹中会产生了许多的**.gnu**文件可以用来绘图,与上面相同的思路,此时递归遍历每个文件夹并执行**plot.sh**这个绘图脚本
```shell
#!/bin/bash 
function getdir(){
    for element in `ls $1`
      do
        dir_or_file=$1"/"$element
    if [ -d $dir_or_file ]
      then
        getdir $dir_or_file
      else  # 下面的全是文件
	  	if [ "${dir_or_file##*.}"x = "gnu"x ];then	# 筛选处特定后缀的文件
    		dir_name=`dirname $dir_or_file` # 读取目录
			file_name=`basename $dir_or_file .gnu` # 读取以.gnu结尾的文件名
			out_file_name="$dir_name/$file_name"  # 定义编号成功的文件名称
			#ifort -mkl $dir_or_file -o $out_file_name.out  # 编译后文件名以out结尾
			dir1=`dirname $out_file_name`
			#echo $dir1
			cd $dir1  # 切换到具体的文件夹
			gnuplot $dir_or_file &
		fi
        #temp_file=`basename $dir_or_file  .f90` #将文件名后缀删除
        #ifort -mkl $dir_or_file -o $temp_file.out  # 编译后文件名以out结尾
        #echo $dir_or_file     # 这里的变量时完整的路径名
    fi
done
}
 
# fold="/home/yxli/te"
fold=`pwd`
getdir $fold
```
最终就可以将所有文件夹下面的程序执行完毕,并进行绘图了.

其实这里还有一个小的缺陷,就是必须等到程序执行完毕之后才会去画图,这当然是必然的事情,程序不执行完毕,没有数据肯定不能绘图,但我们不知道程序何时会执行完毕,所以此时还是需要认为的等待程序执行结束后,才能运行批量绘图的这个脚本,我想可以进行升级,让服务器自动监测程序时候执行完毕,如果执行完毕就可以自动执行绘图的脚本,这样只需要服务器自行判断即可,我们需要的只是去做别的事情,然后再完成之后,去检查计算结果即可,节省时间,高效率才是搞科研的终极目标.
{:.warning}

上面提到的小缺陷我已经在[监测程序运行的Shell脚本](https://yxli8023.github.io/2021/05/15/Shell-Monitor.html)这篇博客中给出了解决方案.
{:.success}

## 小改动
这里稍微修改了一下程序编译后的名称,干脆用数字表示了,防止因为原来的文件名太长,导致可执行文件名也很长,不方便查看.
```shell
#!/bin/bash 
# 递归搜寻文件夹下面所有的.f90或者.f后缀结尾的文件,并利用ifort编译该文件,然后执行
cnt=0  # 定义一个变量,来统计文件夹下面对应程序的个数
function getdir(){
    for element in `ls $1`
      do
		out="out"
        dir_or_file=$1"/"$element
    if [ -d $dir_or_file ]
      then
        getdir $dir_or_file
      else  # 下面的全是文件
		rm *dat *gnu *png 1>/dev/null 2>/dev/null
	  	if [ "${dir_or_file##*.}"x = "f90"x ]||[ "${dir_or_file##*.}"x = "f"x ];then	# 筛选处特定后缀的文件
			cnt=$[$cnt+1]
    		dir_name=`dirname $dir_or_file` # 读取目录
			file_name=`basename $dir_or_file .f90` # 读取以f90结尾的文件名
			#out_file_name="$dir_name/$file_name"  # 定义编号成功的文件名称
			out_file_name="$dir_name/$cnt"  # 定义编号成功的文件名称
			ifort -mkl $dir_or_file -o $out_file_name.$out  # 开始编译fortran文件,编译后文件名以out结尾,以数字命名
			dir1=`dirname $out_file_name`
			#echo $out_file_name.$out
			#echo $dir1
			cd $dir1  # 切换到具体的文件夹,是为了在具体的文件夹中,运行编译好的可执行文件
			#echo $cnt.$out
			./$cnt.$out 1>/dev/null 2>/dev/null &  # 执行该文件夹下面编译好的文件
			#rm $out_file_name.$out
		fi
        #temp_file=`basename $dir_or_file  .f90` #将文件名后缀删除
        #ifort -mkl $dir_or_file -o $temp_file.out  # 编译后文件名以out结尾
        #echo $dir_or_file     # 这里的变量时完整的路径名
    fi
done
}
 
root_dir=`pwd`
getdir $root_dir

```

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)