---
title: wannier90_hr数据的可视化分析
tags: Topology
layout: article
license: true
toc: true
key: a20210310a
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
通常会遇到WannierTools计算紧束缚模型的问题,有时候也会遇到只有wannier90_hr.dat这个数据,那么怎么样可以通过这个数据来分析出对应的紧束缚模型到底是什么样的,这里就利用Mathematica来计算一下如何将wannier90_hr.dat读入,然后可视化的将(x,y,z)方向上的对应的矩阵元素表示出来,这样可以便于问题的分析.
<!--more-->
# Pauli矩阵构建
首先我们先将需要用到的Pauli矩阵构建出来,我这里假定只有两个泡利矩阵的直积,更多自由度之间的直积,在这里进行拓展即可

![png](/assets/images/wannierTools/wannier1.png)


# 数据读入可视化
这里其实就是简单的将数据读入,然后利用表格的方式进行可视化,让整个结构看起来比较清晰,

![png](/assets/images/wannierTools/wannier2.png)

结果如下所示

![png](/assets/images/wannierTools/wannier3.png)

这里的红色标记的正好就是[Wanniertools Tight Binding](https://yxli8023.github.io/2021/03/10/WannierTools-Tight-Binding.html)中对应的(x,y,z),剩下的就是矩阵不分了,接下来的工作就是将这个矩阵表示为Pauli矩阵的直积形式,这个工作在前面我已经涉及到过.

如果想要得到一个对应的矩阵如何分解成Pauli矩阵的直积形式,可以参考[有效边界理论(spinor部分)](https://yxli8023.github.io/2021/01/22/Effective-Edge-Theory-spinor.html)这篇文章.
{:.success}

# 代码下载
所有的内容可以[点击这里下载](/assets/data/wannier90hr.nb)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)