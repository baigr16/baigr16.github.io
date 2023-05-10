---
title: 整理便捷的shell脚本及命令
tags: Code Shell
layout: article
license: true
toc: true
key: a20210517
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
平时在服务器上会用到很多shell命令来进行文件处理等操作,这篇博客主要就是整理我平时学到和用到的方便的shell脚本以及一些方便快捷的命令,方便自己使用的时候进程查阅.
{:.info}
<!--more-->
# 批量kill进程
```shell
#! /bin/bash
com1=$(ps -u yxli|grep out|awk '{print $1}')
#com1=$(ps -u yxli)
for i in $com1
do 
	#echo $i
	kill -9 $i
done
```
从自己(yxli)的进程中(`ps -u yxli`)寻找`out`结尾的进程(`grep out`)并提取第一列的进程号(`awk '{print $1}'`),通过一个`for`循环遍历所有的进程号,然后强制杀死进程(`kill -9 $i`)

这个命令只使用与自己组里的小服务器上使用,如果是服务器集群,或者超算上面,还是乖乖用作业提交脚本和队列系统管理来取消自己提交的错误作业.
{:.success}

# 批量创建文件夹
```shell
for ((i=0;i<=50;i=i+5))do mkdir dir_name$i;done
```
![png](/assets/images/shell/shell1.png)

# 批量复制文件
```shell
for ((i=5;i<=50;i=i+5))do cp bulkimp0.f90 bulkimp$i.f90;done
```
# 批量修改文件内容
```shell
for ((i=0;i<=50;i=i+5))do sed -i "47s/0/$i/g" bulkimp$i.f90;done
```
这个命令会通过循环,将每个文件中第**47行(47s)**的0替换为循环中的**i**.
# 查看用户所有进程
```shell
ps -u yxli
```
这里查询的是用户`yxli`当前所有进程

![png](/assets/images/shell/shell2.png)

```shell
ps|awk '{print $1}'
```
通过通道命令只截取进程ID这一列的内容

![png](/assets/images/shell/shell3.png)

同样可以查询某个用户特定名称的进程,即前面所示的`CMD`这一列中满足特定名称的进程
```shell
ps -u yxli|grep out|awk '{print $1}'
```
这一行可以提取用户`yxli`进程名中包含`out`的所有进程的ID号,如果进程有问题,可以通过前面的批量操作来将所有的进程一次性kill.

# vasp删除多余文件
```shell
rm CHG IBZKPT DOSCAR O* XDATCAR WAVECAR noh* vasprun.xml EIGENVAL REPORT PCDAT PROCAR C*
```
vasp主要的文件只有4个,有时候我们想删除其他的文件,可以用上面的命令,可将其写成脚本,然后再.bashrc中通过一个简单的命令来执行这个操作.

# vasp运行
```shell
alias v5='nohup mpirun -np 10 v5_gam &'
alias v5s='nohup mpirun -np 10 v5_ncl&'
```
开启10个核来执行vasp

# 文件大小显示
```shell
alias ll='ls --block-size=M -l -t'
```
以`M`为单位显示文件大小.

# 批量编译并执行程序
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
将这个脚本文件放在放置在和需要编译的fortran文件相同的文件夹中,这个脚本会自动搜寻当前文件夹中所有后缀名为.f90的文件,然后编译并执行,这里采用的ifort进行编译的,而且因为对角化用到了mkl库中的cheevd这个函数,所有编译选项有-mkl这个参数.
{:,success}

上面的代码可以参考[Fortran + Gnuplot 批量数据输出绘图](https://yxli8023.github.io/2021/05/14/Fortran-Gnu.html)这篇博客的内容.

# 批量执行gnuplot绘图
```shell
#!/bin/sh  
#============ get the file name ===========  
Folder_A=$(pwd) 
for file_a in ${Folder_A}/*.gnu
do 
	gnuplot $file_a  
done
```
当文件夹中有很多gnuplot绘图文件的时候,将前一节的代码进行微改就可以得到批量执行绘图的脚本.

# 递归文件夹编译执行Fortran
```shell
#!/bin/bash 
# 递归搜寻文件夹下面所有的.f90或者.f后缀结尾的文件,并利用ifort编译该文件,然后执行
function getdir(){
    for element in `ls $1`
      do
		    out="out"
        dir_or_file=$1"/"$element
    if [ -d $dir_or_file ]
      then
        getdir $dir_or_file
      else  # 下面的全是文件
		rm *dat *gnui 1>/dev/null 2>/dev/null
	  	if [ "${dir_or_file##*.}"x = "f90"x ]||[ "${dir_or_file##*.}"x = "f"x ];then	# 筛选处特定后缀的文件
    		dir_name=`dirname $dir_or_file` # 读取目录
			file_name=`basename $dir_or_file .f90` # 读取以f90结尾的文件名
			out_file_name="$dir_name/$file_name"  # 定义编号成功的文件名称
			ifort -mkl $dir_or_file -o $out_file_name.$out  # 开始编译fortran文件,编译后文件名以out结尾
			dir1=`dirname $out_file_name`
			#echo $dir1
			cd $dir1  # 切换到具体的文件夹
			./$file_name.$out 1>mes 2>bad &  # 执行该文件夹下面编译好的文件
			# ./$out_file_name.out 1>mes 2>bad &
			rm $out_file_name.$out
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
前面的编译脚本只能在当前文件夹中由Fortran程序的时候运行,如果当前文件夹下面还有文件夹,其中又包含了Fortran程序,这个时候就需要递归的寻找所有的内容,然后进行编译,此时就需要用到上面的脚本.这里可以通过
```shell
 out="out"

 ifort -mkl $dir_or_file -o $out_file_name.$out  # 开始编译fortran文件,编译后文件名以out结尾
```
这两个来边界的修改所有的编译后的可执行程序的名称,我这里采用的是用原来的文件名,将后缀名改为**.out**,即**filename.out**来作为Fortran编译之后的可执行文件名的名称,若想修改,仅以通过修改
```shell
out="out"
```
这一项就可以,也就是通过修改后缀名,可不是去直接修改整个名称,这样可以方便我们知道哪个可执行程序对应的Fortran源代码是哪一个文件.
# 递归文件夹执行gnuplot绘图
```shell
#!/bin/bash 
# 递归搜寻文件夹下面所有的.gnu后缀结尾的文件,并利用gnuplot执行该文件
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
			gnuplot $dir_or_file &  # 执行搜寻到的.gnu后缀的文件
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
这个脚本是在递归编译执行fortran的基础上修改的,只是简单的替换了`gnuplot $dir_or_file & `这个执行命令,以及修改寻找`if [ "${dir_or_file##*.}"x = "gnu"x ];then`后缀名为`.gnu`的文件.

# 批量创建文件并写入内容
```shell
for i in  5.2 5.3 5.4 5.5 5.6 5.7 5.8 ; do
cat >POSCAR$i <<!
cubic diamond
   $i 
 0.0    0.5     0.5
 0.5    0.0     0.5
 0.5    0.5     0.0
  2
Direct
 -0.125 -0.125 -0.125
  0.125  0.125  0.125
!
echo "a = $i"

done
```
这个脚本可以通过循环创建不同的POSCAR文件,而且每个文件中的第二行的内容都不同(`$i`会随着循环改变).代码第二行中的`<<!`表示下一次遇到感叹号`!`就是文件内容的截止位置.

![png](/assets/images/shell/shell-file.png)

# 服务器文件传输
最近需要用命令行来进行一些文件的传输,这里就整理一下如何使用尽量便捷的脚本来实现文件传输
```shell
scp -r name@100.100.100.100:$* /User/home
```
将这个命令写成一个文件`dir-trans.sh`,之后只需要使用
```shell
sh dir-trans.sh file
```
这时候在服务器在执行的时候就会用file替换前面的`scp -r name@100.100.100.100:$* /User/home`中的`*`,这样就可以便捷的传输文件,而不需要每次都需要输入用户名和服务器的ip,同样可以将这个过程再进一步简化,在`.bashrc`中设置一个简单的命令来执行
```shell
alias dir='sh dir-trans.sh'
```
这样的话就只需要执行
```shell
dir file
```
就可以直接通过命令行这种方式来实现和服务器之间的文件传输了.

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)