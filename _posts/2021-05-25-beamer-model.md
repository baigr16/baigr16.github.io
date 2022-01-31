---
title: 华南师范大学Latex幻灯片模板
tags: Study 
layout: article
license: true
toc: true
key: a20210525
pageview: true
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
因为最近正好在写毕业论文,在整理论文的过程中也需要总结一下内容,答辩的时候是需要的,所以就整理了一个符合自己学校风格的`Latex`模板,希望可以在之后的答辩中给自己加点分,下学期的组会报告也可以用这个模板很好的整理一下.
{:.info}
<!--more-->
# 前言
虽然我之前已经在[Latex PPT模板及笔记模板](https://yxli8023.github.io/2020/12/28/note-ppt.html)这篇博客中整理了一个`Latex`写的PPT模板,但使用的是传统了样式,虽然干净整洁,但是显示不出特色,所以就像根据自己学校的校徽配色,来整理一个符合自己学校校徽风格的模板,这样也算是自己对学校做点贡献吧.

我是利用[Charlie Li](https://github.com/CharlieLeee)的[My_Beamer_Template](https://github.com/CharlieLeee/My_Beamer_Template)的模板进行修改的,其实主要就是调整了配色方案和校徽,其余的内容并没有做太多的改动,因为自己为`Latex`写包的代码并不懂,只能根据自己的基本理解和探索进行修改.

# PPT展示
这里利用图片的方式展示一下我制作的PPT模板

![png](/assets/images/latex/beamer_Page1.png)

![png](/assets/images/latex/beamer_Page2.png)

![png](/assets/images/latex/beamer_Page3.png)

![png](/assets/images/latex/beamer_Page4.png)

![png](/assets/images/latex/beamer_Page5.png)

![png](/assets/images/latex/beamer_Page6.png)

![png](/assets/images/latex/beamer_Page7.png)

![png](/assets/images/latex/beamer_Page8.png)

![png](/assets/images/latex/beamer_Page9.png)

![png](/assets/images/latex/beamer_Page10.png)

# 加时钟
在用这个PPT做毕业答辩的时候,因为有时间要求,所以就想在里面增加一个时钟,这样可以让自己在答辩的时候知道时间,不会因为超时每老师嫌弃.
```latex
\usepackage[font=TimesI,timeinterval=1]{tdclock} % 时钟
\date[\initclock\factorclockfont{2.0}\tdtime]{\large 5月25日,2021年} % 显示当前时间
```
增加一个宏包`tdclock`,之后在`\date`处修改一下即可.

**如果使用了这个模板,祝大家答辩顺利.**

# 下载
所有的源代码可以[点击这里下载](/assets/pdf/beamer-model.zip)