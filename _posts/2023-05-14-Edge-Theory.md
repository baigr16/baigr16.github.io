---
title: 拓扑绝缘体(BHZ模型)边界态理论计算
tags:  Topology
layout: article
license: true
toc: true
key: a20230514
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
解析求解拓扑绝缘体边界态
{:.info}
<!--more-->

# 前言

拓扑物态的研究已经成为了凝聚态研究中重要的一部分，我先不赘述与拓扑相关的基本知识了，其实在之前的[量子几何张量与量子度规(1)](https://yxli8023.github.io/2023/05/06/Quantum-metric.html)这篇文章中，已经给出了基本的Berry位相以及Berry曲率的概念。    

2D拓扑绝缘体是受到时间反演对称性保护的具有helical的边界态，在六角晶格上的就是Kane-Mele模型，而在四方晶格上就是BHZ模型来研究。在这里就先关注一下比较简单的BHZ模型，关于数值如何计算边界态在这里就先不关注了，这里主要是通过解析的方式给出拓扑绝缘体的边界态。

![png](/assets/images/20230514/Edge%20Theory_page-0003.jpg)

![png](/assets/images/20230514/Edge%20Theory_page-0004.jpg)

![png](/assets/images/20230514/Edge%20Theory_page-0005.jpg)

![png](/assets/images/20230514/Edge%20Theory_page-0006.jpg)

![png](/assets/images/20230514/Edge%20Theory_page-0007.jpg)


以上就是BHZ模型解析求解得到的边界态，可以看到在能带反转点处，它与数值的结果是一致的。上面给出的是在一个确定了$x，y$方向之后，给出的边界态理论，其实也可以将模型转换到$2D$极坐标系统中，或者令直角坐标轴转动起来，此时$x$和$y$都是角度依赖的，同样可以在这样的情况下给出完整的边界态理论。

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)