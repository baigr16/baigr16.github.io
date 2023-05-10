---
title: VASP+Wannier90编译计算紧束缚能带
tags: Code 
layout: article
license: true
toc: true
key: a20210314
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
最近在学习利用Wannier90结合VASP来做计算,这里先整理一下如何把VASP与Wannier90进行接口,在利用Vaspwiki上的一个例子来测试一下编译好的VASP如何得到紧束缚模型的一些数据.
<!--more-->
这里所有的过程都是参照[一文搞定VASP+wannier90构建紧束缚模型](https://mp.weixin.qq.com/s/bMol75R3qobkbEvMeLQWEg)这篇文章中提供的方法进行的. 这里要强调的事情是, 这个方法只是VASP5.4.4与Wannier90-2.1这个版本进行接口编译的方法, 对其他的版本可能并不适用.
# Wannier90安装
> tar -zxvf wannier90-2.1.0.tar.gz #首先，解压安装包
> cd wannier90-2.1.0/ # 其次，进入文件夹
> cp config/make.inc.ifort make.inc # 然后，准备编译文件(这里老王用的是ifort，注意要检查ifort和mpiifort执行命令)
> make # 接着，编译

完成后，编译库，得到libwannier.a文件

> make lib

这个得到的链接文件,在之后VASP进行编译的时候需要使用,编译完成后如下图所示

![png](/assets/images/vasp/vw1.png)

因为2.1版本是可以与VASP结合的最新版本(注: VASP最新6.2版本已经支持wannier90 3.0版本)。但是2.1版本默认安装与VASP接口并不好，主要是借助肖承诚博士写了一个Fortran的接口(https://github.com/Chengcheng-Xiao/VASP2WAN90_v2_fix)，需要注意的是这个接口是针对VASP 5.4.4版本的。1.2版本与VASP的接口是好的，所以1.2版本默认安装就好。表1对比了两个版本的不同，1.2版本的主要缺点就在于不能构建向上自旋和向下自旋的能带，也就是没有自旋轨道耦合作用的铁磁和反铁磁体系。(至少老王没有想到好办法可以实现)
{.warning}

# VASP编译
首先拷贝VASP2WAN90_v2_fix接口文件中的mlwf.patch 到VASP代码Src目录上一级目录下

![png](/assets/images/vasp/vw2.png)

然后执行如下命令

> patch -p0 < mlwf.patch

接着在VASP makefile.include 文件中加入下面两行,注意路径

![png](/assets/images/vasp/vw3.png)

从路径下面复制编译执行文件,我这里选择的时intel编译器

> cp arch/makefile.include.linux_intel

最后编译即可

> make all

![png](/assets/images/vasp/vw4.png)

这里再bin目录下其实应该是**vasp_std**, 我这里自己修改为了**ab**.

# Si计算实例
这里我是学习了Vaspwiki上的[Bandstructure of Si in GW (VASP2WANNIER90)](https://www.vasp.at/wiki/index.php/Bandstructure_of_Si_in_GW_(VASP2WANNIER90))这个实例, 具体的计算过程可以参考网上的过程, 我这里就只是提供一下我自己的计算结果

![png](/assets/images/vasp/vw5.png)

![png](/assets/images/vasp/vw6.png)

完整的计算过程可以[点击这里下载](/assets/data/Si.zip)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)
