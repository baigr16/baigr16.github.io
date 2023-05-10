---
title: Fortran 并行运算
tags: Code Fortran
layout: article
license: true
toc: true
key: a20210603
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
这篇博客中主要整理一下如何利用Fortran来进行并行运算,这样可以在平时的计算中充分发挥计算机的优势.
{:.info}
<!--more-->
虽然前面也有几个博客整理了python以及Julia的并行计算,但是Fortran的并行计算自己也一直在摸索中,因为这里会涉及到进程之间通讯的问题,正好最近得到了别人写的一段并行代码,通过研究还是很快的搞清楚了如何利用Fortran来进行并行运算,这里就把代码和自己的一些笔记整理出来.

# 并行思想
先来认真的描述一下我在整理代码过程中的,怎么把一个程序分拆成并行.

假设我有一个很多层循环的程序,通常程序执行时间长,也都是遇到重复执行结构,但是整体执行时,总会有一个变量,需要在某一个范围内去执行,假设$x\in[a,b]$,那么我们把这个区间分叉成很多段,$[a,x_1],[x_1,x_2],\cdots,[x_n,b]$,这个时候我们把每一段分别放到一个进程上去执行,最后再把所有进程上执行得到的结果整理到一起,那么我们最终得到的就是在区间$[a,b]$计算的结果了.
{:.success}
这里在进行区间分拆的时候,需要将分拆区间的数目和程序执行时可调用的进程数对应起来,这样才可以保证每个进程上都实行不同区间段的计算,最后把所有的区间段计算得到的结果整理起来才可以正好是完整区间计算的结果.
{:.warning}
这里要对多重循环的结构进行并行,是因为fortran的串行执行的时候,只会调用单个进程来计算,速度上会非常慢,如果可以将计算机上可用进程全部调动起来,就会大大提高执行效率.

# 源代码
```fortran
    program main
    include 'mpif.h'
    !---------------------------------------
    integer  my_rank,p,ierr
    integer i,npc,ntotal
    parameter(ntotal = 30016)
    real*8 re(ntotal)
    real*8 firstRx,delta
    !-----------------------------------
    call MPI_INIT(ierr)     ! 初始化进程
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr) ! 得到本进程在通信空间中的rank值,即在组中的逻辑编号(该 rank值为0到p-1间的整数,相当于进程的ID。)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr) !获得进程个数 size, 这里用变量p保存
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !------------------------
    ! write(*,*)"mu_rank = ",my_rank
    ! write(*,*)"Number of process = ",p
    npc = anint(1.0*ntotal/p)  ! 四舍五入取整
    delta = 0.50d0

    ! 这里通过这个设置,使得不同进程执行不同区间的计算,最后通过MPI_GATHER来将所有的结果收集到一起,得到最终的计算结果
    ! 不同的进程对应的my_rank这个值不同,从而反应在firstRx这个参数上,
    ! firstRx = my_rank*npc*delta + 15.0d0
    firstRx = my_rank*npc*delta
    
    ! call rkky(delta,firstRx,npc,re)  ! 该子过程中,re是返回的结果
    call test(delta,firstRx,npc,re)
    
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    ! 将所有进程结果re(第一个参数)都收集到re(第四个参数)中,此时不同进程中数据长度应该相同
    call MPI_GATHER(re,npc,MPI_real8,re,npc,MPI_real8,0,MPI_COMM_WORLD,ierr)

    !-----------------------------------------
    if (my_rank .EQ. 0) then
    open(1,file='result.dat')
    do i = 1,npc*p
        write(1,*),re(i)
    end do
    endif
    close(1)
    call MPI_FINALIZE(ierr)   ! 结束所有进程
    END PROGRAM 
!======================================================================================
    subroutine test(delta,firstRx,npc,re)
    implicit none
    integer npc,nthe,i1,i2,i3,i4
    real*8 firstRx,re(npc),dthe,delta,R
    real*8 re1
    parameter(nthe=100)
    dthe = 3.1415926d0/nthe
    R = firstRx + delta*npc
    re1 = 0
    ! write(*,111)firstRx,1.0*npc,R
    do i1 = 1,npc
        do i2 = 1,npc
            do i3 = 1,npc
                do i4 = 1,npc
                    re1 = re1 + cos(sin(1.0*i1) + cos(2.0*i2) + sin(1.0*i3) + cos(1.0*i4))
                end do
            end do
        end do
        re(i1) = cos(firstRx)*sin(1.0*i1) + re1
    end do
111 format(3f16.2)
    return
    end subroutine test
!======================================================================================  
    subroutine rkky(delta,firstRx,npc,re)
    implicit none
    integer i,ii,m,n,j,nthe,npc,nk1,nk2,nk3,nk4,nk5,nk6,nx,nqx
    real*8 delta,firstRx,re(npc)
    real*8 R,re11,re1,re2,re3,k0
    real*8 f1,f2,dthe,dk1,dk2,dk3,dk4,dk5,dk6,dx,dqx,qx,x,k,the

    parameter(nqx=7500,nx=2000,nk1=2000,nk2=85,nk3=120,nk4=180,nk5=20,nk6=45,nthe=100)
    parameter(k0=1.0d0)   

    dthe = 3.1415926d0/nthe
    ! 对不同的积分区间选择不同的间隔
    dk1 = 5.0d0/nk1
    dk2 = 25.0d0/nk2
    dk3 = 70.0d0/nk3
    dk4 = 900.0d0/nk4
    dk5 = 1000.0d0/nk5
    dk6 = 8000.0d0/nk6

    dx = 1.0d0/nx
    dqx = 30.0d0/nqx

    do i = 1,npc
        R = firstRx + delta*(i-1)   ! 不同的进程这个R的值是不同的
        re1 = 0.0d0
        do ii = 1,nqx
            qx = dqx/2.0d0 + dqx*(ii-1)
            re11 = 0.0d0

            do m = 1,nx
                x = dx/2.0d0 + dx*(m-1)
                re2 = 0.0d0
            !---------------------------------------------------------
                do n = 1,nk1
                    k = dk1/2.0d0 + dk1*(n-1)
                    re3 = 0.0d0

                    do j = 1,nthe
                        the = dthe/2.0d0 + dthe*(j-1)
                        f1 = (1-x)*(k**2-k0**2)**2 + x*(2*qx*k*cos(the) + qx**2 + k**2-k0**2)**2
                        f2 = qx**2 + 2*qx*k*cos(the)
                        re3 = re3 + (k*f2*(k**2-k0**2)/f1)*cos(qx*R)*exp(-qx**2/64.0d0)*dthe               
                    end do 
                    re2 = re2 + re3*dk1
                end do
            !-------------------------------------------------------------
                do n = 1,nk2
                    k = 5.0d0 + dk2/2.0d0 + dk2*(n-1)
                    re3 = 0.0d0

                    do j = 1,nthe
                        the = dthe/2.0d0 + dthe*(j-1)
                        f1 = (1-x)*(k**2-k0**2)**2 + x*(2*qx*k*cos(the) + qx**2 + k**2-k0**2)**2
                        f2 = qx**2 + 2*qx*k*cos(the)
                        re3 = re3 + (k*f2*(k**2-k0**2)/f1)*cos(qx*R)*exp(-qx**2/64.0d0)*dthe                
                    end do 
                    re2 = re2 + re3*dk2
                end do
            !-----------------------------------------------------------------
                do n = 1,nk3
                    k = 30.0d0 + dk3/2.0d0 + dk3*(n-1)
                    re3 = 0.0d0

                    do j = 1,nthe
                        the = dthe/2.0d0 + dthe*(j-1)
                        f1 = (1-x)*(k**2-k0**2)**2 + x*(2*qx*k*cos(the) + qx**2 + k**2-k0**2)**2
                        f2 = qx**2 + 2*qx*k*cos(the)
                        re3 = re3 + (k*f2*(k**2-k0**2)/f1)*cos(qx*R)*exp(-qx**2/64.0d0)*dthe                
                    end do 
                    re2 = re2 + re3*dk3
                end do
            !------------------------------------------------------------------------
                do n = 1,nk4
                    k = 100.0d0 + dk4/2.0d0 + dk4*(n-1)
                    re3 = 0.0d0

                    do j = 1,nthe
                        the = dthe/2.0d0 + dthe*(j-1)
                        f1 = (1-x)*(k**2-k0**2)**2 + x*(2*qx*k*cos(the) + qx**2 + k**2-k0**2)**2
                        f2 = qx**2 + 2*qx*k*cos(the)
                        re3 = re3 + (k*f2*(k**2-k0**2)/f1)*cos(qx*R)*exp(-qx**2/64.0d0)*dthe                
                    end do 
                    re2 = re2 + re3*dk4
                end do
            !----------------------------------------------------------------------
                do n = 1,nk5
                    k = 1000.0d0 + dk5/2.0d0 + dk5*(n-1)
                    re3 = 0.0d0

                    do j = 1,nthe
                        the = dthe/2.0d0 + dthe*(j-1)
                        f1 = (1-x)*(k**2-k0**2)**2 + x*(2*qx*k*cos(the) + qx**2 + k**2-k0**2)**2
                        f2 = qx**2 + 2*qx*k*cos(the)
                        re3 = re3 + (k*f2*(k**2-k0**2)/f1)*cos(qx*R)*exp(-qx**2/64.0d0)*dthe                
                    end do 
                    re2 = re2 + re3*dk5
                end do
            !-----------------------------------------------------------------------
                do n = 1,nk6
                    k = 2000.0d0 + dk6/2.0d0 + dk6*(n-1)
                    re3 = 0.0d0

                    do j = 1,nthe
                        the = dthe/2.0d0 + dthe*(j-1)
                        f1 = (1-x)*(k**2-k0**2)**2 + x*(2*qx*k*cos(the) + qx**2 + k**2-k0**2)**2
                        f2 = qx**2 + 2*qx*k*cos(the)
                        re3 = re3 + (k*f2*(k**2-k0**2)/f1)*cos(qx*R)*exp(-qx**2/64.0d0)*dthe                
                    end do 
                    re2 = re2 + re3*dk6
                end do
            !--------------------------------------------------------------------------
                re11 = re11 + re2*dx
            end do 
            re1 = re1 + re11*dqx       
        end do
        re(i) = re1
    end do
    return
    end subroutine 
```

# 代码解释
test这个子过程是拿来做测试说明的,可以看到它是一个多重循环结构,最终计算的结果会存储在re中,作为计算返回的结果.
```fortran
	subroutine test(delta,firstRx,npc,re)
    implicit none
    integer npc,nthe,i1,i2,i3,i4
    real*8 firstRx,re(npc),dthe,delta,R
    real*8 re1
    parameter(nthe=100)
    dthe = 3.1415926d0/nthe
    R = firstRx + delta*npc
    re1 = 0
    ! write(*,111)firstRx,1.0*npc,R
    do i1 = 1,npc
        do i2 = 1,npc
            do i3 = 1,npc
                do i4 = 1,npc
                    re1 = re1 + cos(sin(1.0*i1) + cos(2.0*i2) + sin(1.0*i3) + cos(1.0*i4))
                end do
            end do
        end do
        re(i1) = cos(firstRx)*sin(1.0*i1) + re1
    end do
111 format(3f16.2)
    return
    end subroutine test
```
在主程序中,对不同的进程都会得到不同的**my_rank**这个值,根据不同的**my_rank**可以得到不同的**firstRx**,这个变量最后会作为参数输入到**test**这个自过程中,也正是通过**my_rank,firstRx**这两个变量的结合,来实现前面提到的,将整个计算区间分成一系列离散的区间,而且离散区间的数目会和开启的进程数目相同,从而实现在每个进程上分别计算不同离散区间的结果.
```fortran
	call MPI_INIT(ierr)     ! 初始化进程
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr) ! 得到本进程在通信空间中的rank值,即在组中的逻辑编号(该 rank值为0到p-1间的整数,相当于进程的ID。)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr) !获得进程个数 size, 这里用变量p保存
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !------------------------
    ! write(*,*)"mu_rank = ",my_rank
    ! write(*,*)"Number of process = ",p
    npc = anint(1.0*ntotal/p)  ! 四舍五入取整
    delta = 0.50d0

    ! 这里通过这个设置,使得不同进程执行不同区间的计算,最后通过MPI_GATHER来将所有的结果收集到一起,得到最终的计算结果
    ! 不同的进程对应的my_rank这个值不同,从而反应在firstRx这个参数上,
    ! firstRx = my_rank*npc*delta + 15.0d0
    firstRx = my_rank*npc*delta
	call test(delta,firstRx,npc,re)
```
最后,将所有进程计算得到的结果**re**收集到一起,得到最终的计算结果**re**
```fortran
	! 将所有进程结果re(第一个参数)都收集到re(第四个参数)中,此时不同进程中数据长度应该相同
    call MPI_GATHER(re,npc,MPI_real8,re,npc,MPI_real8,0,MPI_COMM_WORLD,ierr)
```
并将得到的结果进行数据保存
```fortran
    if (my_rank .EQ. 0) then
    open(1,file='result.dat')
    do i = 1,npc*p
        write(1,*),re(i)
    end do
    endif
    close(1)
```

这里要注意,既然是要并行,那么你需要并行的那段程序,一定是串行执行的时候需要耗费一定时间的,否则你就没必要并行.
{:.warning}
比如你有一个简单的循环
```fortran
do i4 = 1,npc
	re1 = sin(1.0*i4)
end do
```
这种就完全不必要并行了,就算是并行了也会报错,因为执行时间太短了,进程间的通讯时间可能都会比执行时间要长.

# 编译执行
我是利用mpif90进行编译的
```shell
mpif90 filename.f90 -o f1.out
```
并行执行
```shell
mpirun -np 10 f1.out &
```
此时会开启10个进程来执行这个程序,一般在并行执行程序的时候,最好知道电脑上可用的进程数目是多少,这样可以将所有的进程都调用起来进行计算.


# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)