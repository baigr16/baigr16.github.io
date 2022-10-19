---
title: 笔记模板(Latex)
tags:  Latex
layout: article
license: true
toc: true
key: a20220919
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
在科研中我通常需要整理一些笔记,虽然也想用iPad来整理,但是排版上总是不尽人意,字写得也不好看,所以干脆就整理一个`Latex`模板专门整理学习笔记.这篇Blog就是将自己的这个模板分享出来,主要还是基于[project-logbook](https://github.com/apalha/project-logbook)进行了一些改造.
{:.info}
<!--more-->

这个笔记的模板就是将[科研日志模板(Latex)](https://yxli8023.github.io/2022/09/09/Research-Log.html)这个模板修改了一下，将其中的一些环境块修改了，最后只保留了以`note`环境块，可以通过对其设置不同的颜色框来达到不同的目的，其实这个模板还是很简单的。
# 高亮块
想要将一些文字高亮起来,这里使用了一个`note`环境来实现
```latex
\begin{note}{这里写title}{blue!50}	
\end{note}
```

```latex
\begin{note}{Question}{red}
\end{note}
```
在使用`note`环境的时候，第一个参数`title`就是用来设置高亮标题的,第二个参数则是控制高亮块的颜色。可以将重要的结论或者要注意的一些内容通过设置不同的颜色来进行标记。

![png](/assets/images/latex/note-1.png)

想要设置不同的颜色，可以通过修改颜色`blue!50`产生不同的颜色块,从而来对不同的内容进行高亮标记。下面给一个我整理的笔记示例。

![png](/assets/images/latex/note-2.png)

除此之外，还可以考虑给自己的笔记增加水印，可以在导言区加入下面的命令，这是用来加入一个图片水印，可以自己设计一个专属图片
```latex
\makeatletter
\@namedef{ver@everypage.sty}{9999/99/99}
\makeatother
\usepackage{everypage-1x}
\usepackage[contents=DRAFT, color=red, opacity=0.2]{background}
\backgroundsetup{scale=0.5, angle=45, opacity = 0.1, contents = {\includegraphics[width=\paperwidth, height=\paperwidth, keepaspectratio]{back}}}
```
不过我尝试了一下，加入图片水印会使得产生的PDF有点大，如果是笔记的最终版，可以考虑加图片背景。效果如下图


![png](/assets/images/latex/note-3.png)

在平时整理笔记的时候，可以选择加入文字水印
```latex
\usepackage{draftwatermark}         % 所有页加水印
%\usepackage[firstpage]{draftwatermark} % 只有第一页加水印
\SetWatermarkText{Draft}           % 设置水印内容
%\SetWatermarkText{\includegraphics{fig/texlion.png}}         % 设置水印logo
\SetWatermarkLightness{0.9}             % 设置水印透明度 0-1
\SetWatermarkScale{1}        
```
同样可以在导言区进行设置即可，效果如下图所示

![png](/assets/images/latex/note-4.png)

# 下载
暂时先将模板整理到这里，我觉得我的需求是满足了，感兴趣的可以[点击这里下载](/assets/data/note-model.zip)