---
title: VASP编译安装
tags: vasp
layout: article
license: true
toc: true
key: a202008092
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
pageview: true
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
假期空闲，又不能外出学习交流，正好趁这段时间入门一下第一性原理计算，首先从安装VASP开始，这里记录一下自己的安装过程，说不定之后还会用的到。
{:.success}
<!--more-->
# intel fortran安装
首先需要安装Intel Parallel studio XE2019，如果有教育邮箱，可以申请试用，有效期是一年，没有的话请自行解决，百度一堆教程。
进入软件文件夹之后会看到一个**install.sh**的文件,执行命令`./install.sh`这里默认你已经拥有root权限。按照软件的提示一路往下走，之后就是要输入序列号了，如果是教育邮箱申请的，那么你的邮箱会收到这个序列号，我的是S477-LJST6J4M。之后的安装同样遵循默认选项，虽然可以自行选择哪些需要安装，哪些不需要安装，但因为这是新手初学，暂时先不关心这个问题。接下来就是一段时间的等待之后，就可以成功安装了。

![png](/assets/images/vasp/intel-1.png)

![png](/assets/images/vasp/intel-2.png)

![png](/assets/images/vasp/intel-3.png)

安装结束之后，并不代表软件可以使用了，如果熟悉Linux系统的话，就知道我们需要设定一下环境变量，才能保证每次进入终端之后可以成功使用intel fortran。
首先需要找到正确的安装路径，我的安装路径为

```shell
/opt/intel
```

找到里面的这个文件 **psxevars.sh**
对应的路径为 /opt/intel/parallel_studio_xe_2019
打开`.bashrc`文件(注意前面有个英文句号，说明这是Linux中的隐藏文件)，然后将环境变量加入

```shell
source /opt/intel/parallel_studio_xe_2019/psxevars.sh
```

![png](/assets/images/vasp/intel-5.png)

加入后保存.bashrc文件，然后执行`source .bashrc`这样就成功的将intel fortran安装到了机器上

![png](/assets/images/vasp/intel-6.png)
# 检查intel fortran是否可以成功安装

```shell
icc -v    ifort -v    这两个命令会分别返回c和fortran编译器的版本
which ifort  这个命令则会返回ifort这个执行命令的路径，如果成功安装，则上面的命令都不会报错
echo $MKLROOT  这个命令是告诉你你的mkl函数库在哪里，这个库函数主要是用来做矩阵运算的，一定要正确安装
```

# 编译intel fftw3
在成功安装好intel之后，接下来就是要在intel安装位置的根目录中找到fftw3，如果是按照上面的方式安装的，那么对应的路径位置是

```shell
/opt/intel/compliers_and_libraies/linux/mkl/interfaces/fftw3xf
```

![png](/assets/images/vasp/intel-4.png)

进入这个目录之后，执行

```shell
make libintel64
```

等待完成之后即可

![png](/assets/images/vasp/intel-7.png){:width="330px",:height="495px"}![png](/assets/images/vasp/intel-8.png){:width="330px",:height="495px"}

# VASP编译
这里先说几句，我也是初学者，但是就我所知VASP是有许多不同的版本的，不同的版本可能编译方法都会有所不同，这可能也是网上有很多教程来教你怎么安装VASP，但是由于你手头不论是VASP或者intel fortran的版本和作者的版本都是不太相同的，所以即使是完整参照了教程，最后发现安装还是有一堆问题出现。我先在用的VASP的版本是5.4.4，Fortran编译器是ifort version 19.0.4.243(ifort -v查询结果)。

将VASP5.4.4解压之后，文件夹结构如下图所示

![png](/assets/images/vasp/vasp-1.png){:width="330px",:height="495px"}![png](/assets/images/vasp/vasp-2.png){:width="330px",:height="495px"}

在arch中有对于不同机器，已经写好的makefile文件，接下来我们要做的就是选择合适自己机器的makefile文件，然后进行编译。这里我们选择**makefile.include.linux_intel**，将这个文件复制到和arch文件夹同级的目录下，如下图所示，然后将这个复制过来的文件重命名为**makefile.include**。接下来就是要对这个文件内容进行一些小小的调整，源文件第20行的内容为

```shell
OFLAG = -02
```

![png](/assets/images/vasp/vasp-3.png)

这里将它修改为
```shell
OFLAG = -02 -xhost
```

增添-xhost参数可以使ifort编译出的程序能够利用当前机子CPU支持的最高档次SIMD指令集。修改之后进行保存，然后就可以开始编译了。

```shell
make all
```

剩下的就是长时间的等待，附一张无意义的图来说明等待是有多么无聊

![png](/assets/images/vasp/vasp-4.png)

# 编译结束
编译完成之后，在你的之前进行操作的目下面，将会出现一个**bin**文件夹，里面会有3个文件，如图所示

接下来将编译好的标准版本的VASP的可执行程序改名为vasp，接着把这个VASP的执行文件加入到环境变量中，在**.bashrc**中加入下面这条命令
>export PATH=$PATH:/sob/vasp.5.4.4/bin

这个命令在不同人的机器上是需要进行修改的，比如**/sob/vasp.5.4.4/bin**，这是我自己机器上编译成功的VASP的可执行文件，这里换成你自己的即可，完成之后记得重新加载一下**.bashrc**，即执行`source .bashrc`，安装到此结束。

# 运行
> mpirun -np 4 vasp(四核并行)

这个命令是需要在一个文件夹下面执行的，文件夹下面需要必备VASP运行的四个文件，这不是这里要讨论的问题。

# 后话
据我了解到，VASP应该还可以去集成别的东西，但是我自己也并不是很了解，所以暂时先把这个留作一个未完成的任务。

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)