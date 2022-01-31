---
title: 监测程序运行的Shell脚本
tags: Code Shell
layout: article
license: true
toc: true
key: a20210515
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
这篇博客整理一个在Linux服务器上通过进程名来监测程序是否执行的脚本,然后根据程序是否执行来进行其他的一些操作,这对在服务器上跑程序,然后自动化执行某些操作还是非常方便的.
{:.info}
<!--more-->
在前面[Fortran + Gnu 批量计算](https://yxli8023.github.io/2021/05/14/Fortran-Gnu.html)这篇博客中,虽然可以批量的进行fortran程序的编译和执行,在程序执行完成之后再批量的进行gnuplot绘图,这里的缺陷就是必须人为的取查看前面编译的fortran程序是否已经执行完毕,谁也不会准确的指导自己的程序何时会在服务器上执行完毕,那么就可以通过一个脚本让服务器自行监测程序的运行,当fortran程序执行完毕之后,接下来就可以让服务器再自行启动gnuplot的绘图程序.

# 程序监测
这篇博客的脚本其实也就是对[Fortran + Gnu 批量计算](https://yxli8023.github.io/2021/05/14/Fortran-Gnu.html)这篇博客的一个缺陷补充,但是监测程序修改一下同样可以监测别的东西,所以就单独整理了出来.
- 初级版
```shell
#!/bin/bash
while :
do 
	sleep  10m # 每10分钟进行一次监测
	COUNT=$(ps -ef |grep "out" |grep -v "grep" |wc -l) # 从进程中查看out结尾的程序数目
	echo $COUNT>>mes.dat  # 打印进程中out结尾的数量
	if [ $COUNT -eq 0 ];then
		echo Not running  >>  mes.dat # 如果out结尾进程数量不为零,则执行这一项
		sh auto-plot-fold.sh &
		exit 0
	else
		echo Is Run >>  mes.dat   # 如果out结尾进程数量为零,则执行这一项	
		#ps -ef|grep "monitor.sh">>mes.dat	
		echo "ID:">>mes.dat
		ps -ef | grep "out" | grep -v grep | awk '{print "ID:" $2 "   Time:" $7}'>>mes.dat
		# exit 0
	fi
done
```

- 完善版

```shell
#!/bin/bash
rm mes.dat
while :
do 
	conut=$(ps -ef |grep "out" |grep -v "grep" |wc -l) # 从进程中查看out结尾的程序数目
	echo ID_Number: $conut>>mes.dat  # 打印进程中out结尾的数量
	current_time="`date +"%Y-%m-%d %H-%M-%S"`" #获取系统当前时间
	if [ $conut -eq 0 ];then
		echo current_time:$current_time >> mes.dat
		#echo Not running  >>  mes.dat # 如果out结尾进程数量不为零,则执行这一项
		sh auto-plot-fold.sh &
		exit 0
	else
		echo current_time:$current_time >> mes.dat
		echo $conut processing is running >>  mes.dat   # 如果out结尾进程数量为零,则执行这一项	
		#ps -ef|grep "monitor.sh">>mes.dat	
		ps -ef | grep "out" | grep -v grep | awk '{print "ID:" $2 "   Time used:" $7}'>>mes.dat
		# exit 0
		echo -e "\n" >> mes.dat
	fi
	sleep  10m # 每10分钟进行一次监测
done
```

## 代码解释
这里用了一个无限循环,每10分钟进行一次监测
```shell
sleep  10m # 每10分钟进行一次监测
```
接下来进行进程识别
```shell
COUNT=$(ps -ef |grep "out" |grep -v "grep" |wc -l) # 从进程中查看out结尾的程序数目
```
这一行代码主要是在进程中**ps -ef**提取以**out**结尾的进程**grep "out"**,然后统计了这种名称的进程的数目**wc -l**,这里面会涉及到一些基本的shell的语法,看不懂请自行补充知识,我就不去介绍了.

接下来就进行判断,因为我在对fortran程序批量编译执行的时候,就已经将所有的执行文件都已**out**结尾,所以上面在识别进程的时候才可以只关注**out**结尾的进程
```shell
if [ $COUNT -eq 0 ];then
```
如果统计的进程数量等于0,说明此时所有的程序都已经执行完毕,就可以让服务器进行绘图程序启动了
```shell
echo Not running  >>  mes.dat # 将一些信息追加到mes.dat这个文件中
sh auto-plot-fold.sh &  # 启动绘图程序
exit 0 # 退出无限循环
```
在启动了绘图程序之后,就不必在进行监测了,可以退出循环.

如果这种进程的数量不为0,那么此时所有的程序还没有执行完毕,那么向**mes.dat**中写入一些信息,方便我们去查看
```shell
echo Is Run >>  mes.dat   # 如果out结尾进程数量为零,则执行这一项	
#ps -ef|grep "monitor.sh">>mes.dat	
echo "ID:">>mes.dat
ps -ef | grep "out" | grep -v grep | awk '{print "ID:" $2 "   Time:" $7}'>>mes.dat
# exit 0
```

以上就是整个代码的解释.

这里有一点需要注意,因为是在执行shell脚本,此时同样可以让这个脚本在后台自动执行,假设这个脚本名为**monitor.sh**,那么需要利用
`nohup sh monitor.sh &`来提交任务,这样即使我们退出了终端,那么这个任务也会在后台自动执行.
{:.warning}

# 升级版
```shell
#!/bin/sh  
#============ get the file name ===========  
rm *.gnu *.png no* *.dat *.out 1>/dev/null 2>/dev/null mes.dat &
Folder_A=$(pwd) 
for file_a in ${Folder_A}/*.f90
do 
	temp_file=`basename $file_a  .f90` 
	ifort -mkl -O3 -CB $file_a -o $temp_file.out 
    nohup ./$temp_file.out 1>/dev/null 2>/dev/null &
done

while :
do 
	sleep  10 # 每10秒钟进行一次监测
	conut=$(ps -ef |grep "out" |grep -v "grep" |wc -l) # 从进程中查看out结尾的程序数目
	echo ID_Number: $conut>>mes.dat  # 打印进程中out结尾的数量
	current_time="`date +"%Y-%m-%d %H-%M-%S"`"
	if [ $conut -eq 0 ];then
		echo current_time:$current_time >> mes.dat
		#echo Not running  >>  mes.dat # 如果out结尾进程数量不为零,则执行这一项
		# 程序执行完毕开始绘图
		Folder_A=$(pwd) 
		for file_a in ${Folder_A}/*.gnu
		do 
			gnuplot $file_a  1>/dev/null 2>/dev/null &
		done
		exit 0 # 绘图结束后退出死循环
	else
		echo current_time:$current_time >> mes.dat
		echo $conut processing is running >>  mes.dat   # 如果out结尾进程数量为零,则执行这一项	
		#ps -ef|grep "monitor.sh">>mes.dat	
		ps -ef | grep "out" | grep -v grep | awk '{print "ID:" $2 "   Time used:" $7}'>>mes.dat
		# exit 0
		echo -e "\n" >> mes.dat
	fi
done
```
这里把绘图与程序编译执行放在了一起,如果在编写fortran程序的时候,同时对相应的数据写出对应的gnuplot绘图代码,那么可以利用这个程序将程序编译执行和作图一起完成,这样可以节省很多时间.

