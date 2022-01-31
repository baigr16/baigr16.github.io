---
title: Mathematica绘制漂亮的散点密度图  
tags: Plot Mathematica
layout: article
license: true
toc: true
key: a20210704
pageview: true
cover: /assets/images/Mma/den2.png
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
在通过计算得到二维密度图之后,如何绘制一幅漂亮的密度图,可以美观的展示自己的结果,这里就来整理一下如何利用Mathematica来绘制这样的密度图.
{:.info}
<!--more-->
- 首先把代码贴出来

![png](/assets/images/Mma/den1-code.png)

- 首先根据颜色进行区分
![png](/assets/images/Mma/den1.png)

- 在颜色区分的基础上,再根据格点的大小来区分密度
![png](/assets/images/Mma/den2.png)

# 3D
- 首先把代码贴出来

![png](/assets/images/Mma/den3-code.png)

- 首先根据颜色进行区分
![png](/assets/images/Mma/den3.png)

- 在颜色区分的基础上,再根据格点的大小来区分密度
![png](/assets/images/Mma/den4.png)

# 代码下载
源代码可以[点击这里下载](/assets/data/mma-density.nb)


