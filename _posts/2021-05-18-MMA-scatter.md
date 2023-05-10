---
title: Mathamatica 绘制旋转的散点密度图
tags: Plot Mathematica
layout: article
license: true
toc: true
key: a20210518
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
cover: /assets/images/Mma/r1.png
author: YuXuan
show_author_profile: true
---
这里整理一下如何利用散点数据绘制密度分布图,由于直接利用Mathematica的`ListDensityPlot`绘制得到的图实在不忍直视,所以改进了一下利用`Graphics`来对数据做一个漂亮的展示,并加上colorbar.
{:.info}
<!--more-->
# 散点密度图
当利用程序求解得到密度图数据后,利用Mathamatica直接绘制得到的效果并不好,这里先争议如何通过`Graphics`来绘制密度图,可以得到很好的显示效果

![png](/assets/images/Mma/den-scatter1.png)

![png](/assets/images/Mma/r1.png)

## 函数封装

![png](/assets/images/Mma/den-scatter2.png)

![png](/assets/images/Mma/r2.png)

这里在绘制的时候,设置了一下旋转,将整个图形转动了45度,这个转动角度时可以调节的,可以通过**RotationMatrix**中的转动角度来进行调整,但是在调整了旋转角度之后会使得colorbar和图形重叠到一起,这个时候需要重新设置**Inset**中的第二个参数,来调整插入的corlorbar的位置,同时也可以调整第四个参数来设置colorbar的尺寸大小.至于配色也可以自行设置,具体颜色可以查阅帮助文档.
{:.warning}

# 代码下载

上面所有的程序以及用到的数据可以[点击这里下载](/assets/data/rotate-den.zip).

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)