---
title: MPI+Fortran 并行编程笔记
tags: MPI Fortran 
layout: article
license: true
toc: true
key: a202106010
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

这篇博客整理一下自己使用MPI对Fortran进行并行化的时候记录的笔记。

{:.info}
<!--more-->

# 前言

之前虽然也整理过Fortran并行的一些东西，但是没有对MPI和Fortran并行有一个很好的理解，这里就整理一下自己在对Fortran程序进行并行化的时候，并顺便学习MPI时候的一些笔记。

# MPI

MPI其实就是一个管理通信结构的一个语言，在很多地方也都有使用到，我自己现在了解到的就是和Fortran还有C语言几何，来编写并行的程序。而且我自己经常使用的Python还有Julia其实也都有相关的MPI通信的库，所以这个还是很有用的，尤其是在计算有些东西的时候，速度太慢就使得每周的组会内容很少，因为算不出来结果来让自己看起来一周都在摸鱼。费话不多少，这里先来简单的介绍几个MPI最基本的函数。

## MPI_INIT() :

这个函数就是用来出事化MPI环境的，它有一个参数返回值，通过这个返回值就可以判断出是否初始化好了MPI的环境，毕竟只有初始化完成了，才能开始进行后面的过程。

```fortran
PROGRAM hello_world_mpi
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror

call MPI_INIT(ierror) 
```

## MPI_COMM_SIZE

说到并行，那么必然就要提到你是多少个核在并行，那么程序如何知道你一共分配了多少个核来执行计算呢，就是通过`MPI_COMM_SIZE`这个函数，我们就可以在程序中获取到我们所提供的核数。作为并行化的程序，自然是你给我多少个核我就开多少个线程，所以这里其实也就是让程序知道一共可以有多少个线程。

## MPI_COMM_RANK

知道了一共有多少个线程之后，也就是知道了可以将自己的任务分裂成多少份，此时就需要将每一份分配到一个线程上进行计算。但是我们还是需要知道到底哪一份任务分类到了哪一个线程上。也就是说每一个线程都会有自己的一个线程`id`，通过这个`id`来判断到底是哪一份程序在哪一个线程上执行。如下面的程序所示

```fortran
PROGRAM hello_world_mpi
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

print *, 'Hello World from process: ', process_Rank, 'of ', size_Of_Cluster
```

这里

```fortran
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
```

其中`size_Of_Cluster`就返回我们提供的计算核数量，

```fortran
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)
```

`process_Rank`则给出每一个进程的`id`编号。有了上面的这几个函数，还需要另外一个结束MPI并行的函数。

## MPI_FINALIZE()

这个函数的作用就是告诉程序，到这里MPI通信就结束了，草率一点说就是结束了Fortran的并行，所以上面的程序完整版应该是

```fortran
PROGRAM hello_world_mpi
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror, tag

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

print *, 'Hello World from process: ', process_Rank, 'of ', size_Of_Cluster

call MPI_FINALIZE(ierror)
END PROGRAM
```

编译并执行

```fortran
mpiifort hello.f90 -o a.out
mpirun -np 4 ./a.out
```

这里就给了4个核进行并行，结果如图所示

![png](/assets/images/Fortran/mpi-1.png)

上面的过成中看到了，我们开启了4个核进行计算，而且从我们前面的讨论并行的时候，也说到了就是将要执行的任务，根据我们所提供的核数分割成不同的部分，每一部分之间都是独立执行的，从上面的执行结果中可以看到，确实如此，因为每个任务都是独立进行的，所以它们的输出顺序是随机的。但是我们通常会遇到一个问题，那就是我们需要这些任务是顺序执行的，也就是`id`靠前的就先执行，让执行过程是循环的，科学用于角同步执行，所以上面的执行其实是异步执行，这里就来看看如何进行同步执行。

{:.success}



## MPI_BARRIER

`MPI_BARRIER`这个函数的功能就是来实现我们上面所说的，让每一份分开执行的程序，是按照一定顺序进行的。将上面的代码进行修改，加入这个函数

```fortran
PROGRAM hello_world_mpi
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

DO i = 0, 3, 1
    IF(i == process_Rank) THEN
        print *, 'Hello World from process: ', process_Rank, 'of ', size_Of_Cluster
    END IF
    call MPI_BARRIER( MPI_COMM_WORLD, i_error)
END DO

call MPI_FINALIZE(ierror)
END PROGRAM
```

这里加入了

```fortran
call MPI_BARRIER( MPI_COMM_WORLD, i_error)
```

这个函数，也就实现了同步执行的目标，结果如图所示

![png](/assets/images/Fortran/mpi-2.png)

# 信息通讯

MPI的主要作用其实就是用来控制各线程之间的通讯，也是它最主要的应用，这里就来看看如何使用MPI的通讯功能，其实主要就是两个函数

```fortran
call MPI_SEND(integer message, integer count, MPI_Datatype datatype, integer dest,
integer tag, MPI_Comm comm, integer ierror);
call MPI_RECV(integer data, integer count, MPI_Datatype datatype, integer from,
integer tag, MPI_Comm comm, MPI_Status* status, integer ierror);
```

这里一个变量一个变量的来分析这两个函数都需要什么样的输入

## MPI_SEND

```fortran
integer message         !第一个参量是我们想要发送的数据
integer count           !第二个参量是我们想要发送数据的数量
MPI_Datatype datatype   !The MPI-specific data type being passed through the array.
integer dest            !第四个参数用来说明我们要向哪一个线程发送数据
integer tag             !这是一个消息标签
MPI_Comm comm           !The MPI Communicator handle.
integer ierror          !An error handling variable.
```

## **MPI_RECV**

```fortran
integer message:        !第一个是我们要接受的数据
integer count:          !第二个参量是要接受数据的个数
MPI_Datatype datatype:  !The MP-specific data type being passed through the array.
integer from:           !第四个参数用来明确接受的数据是从哪一个线程发送过来的
integer tag:            !这是一个消息标签
MPI_Comm comm:          !The MPI Communicator handle.
MPI_Status* status:     !Status object.
integer ierror          !An error handling variable.
```

这里给出一个例子，从一个线程发送数据到另外一个线程

```fortran
PROGRAM send_recv_mpi
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror, message_Item

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

IF(process_Rank == 0) THEN
    message_Item = 42
    call MPI_SEND(message_Item, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, ierror)
    print *, "Sending message containing: ", message_Item
ELSE IF(process_Rank == 1) THEN
    call MPI_RECV(message_Item, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
    print *, "Received message containing: ", message_Item
END IF

call MPI_FINALIZE(ierror)
END PROGRAM
```

这里再来看一下两个通讯函数的具体参数

```fortran
MPI_SEND(
        message_Item,       !存储要发送的数据
        1,                  !明确要发送的数据的个数
        MPI_INT,            !MPI_TYPE of the message we are sending.
        1,                  !这个数据要发送到1这个线程
        1,                  !Message Tag
        MPI_COMM_WORLD      !MPI Communicator
        ierror              !Error Handling Variable
)
MPI_RECV(
        message_Item,       !用来存储接受到的数据
        1,                  !要接受数据的数量
        MPI_INT,            !MPI_TYPE of the message we are sending.
        0,                  !数据是从哪一个线程发送过来的
        1,                  !Message Tag
        MPI_COMM_WORLD      !MPI Communicator
        MPI_STATUS_IGNORE   !MPI Status Object
        ierror              !Error Handling Variable
)
```

代码执行结果为

```shell
Sending message containing: 42
Received message containing: 42
```

# 数据分发与收集

在MPI中有时候利用群组分发也是非常有用的，其实也就会说将一组数据分发到不同的线程中去，或者从不同的线程中将数据进行汇总。

## **MPI_Scatter**

```fortran
integer send_Var            !要分发的数据
integer send_Count          !分发数据的数量
MPI_Datatype send_Type      !MPI Datatype of the data that is scattered.
integer recv_Var            !Variable that will store the scattered data.
integer recv_Count          !每个线程要结束的数据的个数
MPI_Datatype recv_Type      !MPI Datatype of the data that will be received.
integer root_Process        !执行数据分发的线程号
MPI_Comm comm               !The MPI_Communicator.
integer ierror              !An error handling variable.
```

## **MPI_Gather**

```fortran
integer send_Var            !Variable storing the value that will be sent.
integer send_Count          !接受数据的数量
MPI_Datatype send_Type      !MPI Datatype of the data that is sent.
integer recv_Var            !Variable that will store the gathered data.
integer recv_Count          !从每个线程中接受数据的数量
MPI_Datatype recv_Type      !MPI Datatype of the data that will be received.
integer root_Process        !收集数据所对应的线程号，也就是让那个线程来收集数据
MPI_Comm comm               !The MPI_Communicator.
integer ierror              !An error handling variable.
```

这里同样给一个例子来说明这个函数的用法

```fortran
PROGRAM scatter_mpi
include 'mpif.h'

integer process_Rank, size_Of_Cluster, ierror, message_Item
integer scattered_Data
integer, dimension(4) :: distro_Array
distro_Array = (/39, 72, 129, 42/)

call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)
call MPI_Scatter(distro_Array, 1, MPI_INT, scattered_Data, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror);

print *, "Process ", process_Rank, "received: ", scattered_Data
call MPI_FINALIZE(ierror)
END PROGRAM
```

执行结果如下

```shell
Process 1 received: 39
Process 0 received: 72
Process 3 received: 129
Process 2 received: 42
```

# 参考

