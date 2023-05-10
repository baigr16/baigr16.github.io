---
title: 空间群学习工具SpaceGroupIrep
tags: Group-Theory Mathematica
layout: article
license: true
toc: true
key: a20210419
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
pageview: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
这里整理一下学习空间群时候的一个Mathematica软件包[SpaceGroupIrep](https://github.com/goodluck1982/SpaceGroupIrep),在学习THE MATHEMATICAL THEORY OF SYMMETRY IN SOLIDS这本书中的内容的时候, 结合这个工具可以对书中的内容有一个更深入的理解. 这个软件包也是基于这本书中的内容编写的, 而且可以对书中的一些表进行快速查询, 结合起来对学习空间群相关的知识大有帮助.
{:.info}
<!--more-->
# 前言
自己最近一段时间在学习THE MATHEMATICAL THEORY OF SYMMETRY IN SOLIDS(Representation theory for point groups and space groups)这本书, 虽然认真的学习了一遍, 但是对第四章中的抽象表示还是有一些不明白的地方, 正好可以结合SpaceGroupIrep来对书中的内容再进行一遍学习, 相信结合这个软件包, 可以对书中的内容有进一步的理解.

首先将软件包下载, 然后复制到
```shell
$UserBaseDirectory/Applications
```
这个目录中就可以, 记得将下载下来的软件包的名称修改为SpaceGroupIrep, 这样想要使用的话在Mathematica中执行
```shell
<< "SpaceGroupIrep`"
```
就可以.
# 布里渊区(BZ)展示
首先就是展示不同的布里渊区

![png](/assets/images/20210419/B4.png)

![png](/assets/images/20210419/B3.png)

![png](/assets/images/20210419/B2.png)

![png](/assets/images/20210419/B1.png)

# 转动操作矩阵
这里可以直接通过查表的方式来确定不同BZ对应的转动操作, 也就是书中的Table3.2, 利用软件包可以得到这些操作对应的矩阵形式

![png](/assets/images/20210419/B5.png)

# 双群转动操作
同时也可以计算对应双群的操作表示, 对应书中的Table6.1

![png](/assets/images/20210419/B6.png)

# 小群,Herring小群,高对称点群元
通常我们会需要确定空间群高对称点上的操作群元有哪些, 可以利用**getLGElem**来计算. 而且在计算高对称点的不可约表示, 以及空间群的不可约表示的时候, 会用到小群,Herring小群的群元, 这些信息可以通过**getLGElem,getHLGElem,getCentExt**. 这与这些具体是什么意思, 在这个博客中就不解释了, 这里只是简单的记录并整理一下这个软件包工具性的一面.

![png](/assets/images/20210419/B7.png)

# 群元乘法
群元之间的运算关系同样可以通过计算得到, 主要会用到下面几个函数

![png](/assets/images/20210419/B8.png)

# 抽象群
在查表的时候, 抽象群是经常要用到的, 这里可以通过函数来直接查询书中的抽象群

![png](/assets/images/20210419/B9.png)

# 空间群,小群不可约表示

![png](/assets/images/20210419/B10.png)

![png](/assets/images/20210419/B11.png)

![png](/assets/images/20210419/B12.png)

# 空间群直积分解

![png](/assets/images/20210419/B13.png)

# 能带对应小群的不可约表示

在材料能带计算的过程中, 对每条能带标记其不可约表示是比较重要的一个事情, 这个内容利用SpaceGroupIrep同样可以, 只不过需要结合其他的工具, 这部分内容我并不会, 所以这里就只是介绍一下, 具体我也还在学习过程中.

# 代码
所有上面的内容, 可以[点击这里下载](/assets/data/SpaceGroup.zip)



# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)

