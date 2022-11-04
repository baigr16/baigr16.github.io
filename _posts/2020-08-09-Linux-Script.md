---
title: Linux中批量执行编译并运行Fortran
tags: Linux Fortran Shell
layout: article
license: true
toc: true
pageview: true
key: a20200809
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
最近因为要大量重复的跑一些程序，而且只是参数的小修，所以干脆花点时间整理一个界几个shell脚本，来自动的完成程序的编译及执行。
{:.info}
<!--more-->
# 批量编译Fortran并运行
```shell
#!/bin/sh  
#============ get the file name ===========  
Folder="/home/yxli/te"  	#要批量编译哪个文件夹下面的Fortran
for file_name in ${Folder}/*.f90
do 
	temp_file=`basename $file_name  .f90` 
	ifort -mkl $file_name -o $temp_file.out 
	./$temp_file.out &   # 编译成功之后自动运行
done
rm *out   # 删除编译后文件
```
# 递归的读取指定文件夹下面的所有Fortran文件并编译执行
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
	  	if [ "${dir_or_file##*.}"x = "f90"x ]||[ "${dir_or_file##*.}"x = "f"x ];then	# 筛选处特定后缀的文件
    		dir_name=`dirname $dir_or_file` # 读取目录
			file_name=`basename $dir_or_file .f90` # 读取以f90结尾的文件名
			out_file_name="$dir_name/$file_name"  # 定义编号成功的文件名称
			ifort -mkl $dir_or_file -o $out_file_name.out  # 编译后文件名以out结尾
			dir1=`dirname $out_file_name`  # 读取编译成功文件的路径,只提取目录
			cd $dir1  # 切换到具体的文件夹
			./$file_name.out 1>mes 2>bad &  # 执行该文件夹下面编译好的文件
			# ./$out_file_name.out 1>mes 2>bad &
			# rm $out_file_name.out
		fi
    fi
done
}
 
root_dir="/home/yxli/te"
getdir $root_dir
```
# 一些自己定义的Linux命令
```shell
alias ps='ps -u yxli'   # 查看自己所有进程
alias proid='ps|grep *.out'  #查看以out结尾的进程(正好可以配合前面批量编译Fortran来使用)
```