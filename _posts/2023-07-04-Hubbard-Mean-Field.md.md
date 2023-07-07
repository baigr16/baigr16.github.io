---
title: Hubbard 模型平均场(简单介绍)
tags:  Superconductor
layout: article
license: true
toc: true
key: a20230704
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
简单先整理一下Hubbard模型平均场处理的方法。
{:.info}
<!--more-->

# 前言
Hubbard模型作为凝聚态研究中许多有趣效应的出发点，具有很重要的研究价值。在固体物理最简单的处理中，就是忽略了电子与电子之间的相互作用，也就有了自由电子近似等处理方法。因为在不考虑电子之间相互作用的时候，哈密顿量是两算符形式，可以将其表示为矩阵形式，从而严格求解。但电子之间的相互作用总是存在的，有时候不可避免的是要考虑相互作用的，此时Hubbard相互作用会引入四算符形式，哈密顿量不能再向两算符形式那样简单的表示为矩阵，进行对角化研究该模型。研究超导必然会接触Hubbard模型，有时候则是通过平均场的方法进行处理，该方法简单(目前我也只会这个，强关联太难了，暂未涉猎)，本质上的想法就是将原本的四算符形式分解成两算符再乘以另外两个算符的平均值。从处理方法上来看还是比较简单的，但也的确在研究一些问题的时候很好的抓住了问题的核心，比如用平均场方法研究超导时，将Hubbard相互作用分解成两个产生算符乘以两个湮灭算符的平均值，这样就将哈密顿量变成了都是有两个算符组成的形式(在超导里面这样得到的就是BdG哈密顿量，可参考正常态哈密顿量到BdG哈密顿量的构建这篇文章)。在研究比如拓扑超导的时候，其实出发点也都是BdG哈密顿量，只不过其与绝缘体哈密顿量相比较而言，超导打开的能隙就类比于绝缘体的能隙，如果只从形式上来，似乎也没什么区别。比如看拓扑绝缘体与2D拓扑p波超导体，看起来二者的哈密顿量几乎是一样的，只不过就是其中Pauli矩阵的自由度代表了不同的含义，在超导中通常就说是Nambu自由度(粒子-空穴)。而且利用平均场研究研究其它的序也是可以的(铁磁序，自旋密度等)。这里先简单的介绍一下Hubbard模型平均场处理的内容，关于Hubbard模型更加详细的可以参考The Hubbard Model: A Computational Perspective和The Hubbard Model这两篇在Annual Review of Condensed Matter Physics上的综述文章。
{:.info}


![png](/assets/images/20230704/model_page-0003.jpg)

![png](/assets/images/20230704/model_page-0004.jpg)

![png](/assets/images/20230704/model_page-0005.jpg)

![png](/assets/images/20230704/model_page-0006.jpg)



# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)