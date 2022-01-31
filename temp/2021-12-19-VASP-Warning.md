---
title: VASP报错及修复
tags: Topology 
layout: article
license: true
toc: true
key: a20211219
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
这篇博客主要整理一下自己在VASP使用和学习过程中遇到的错误以及如何修复这些报错.
{:.info}
<!--more-->
> VERY BAD NEWS! internal error in subroutine IBZKPT:
> Reciprocal lattice and k-lattice belong to different class of lattices. Often results are still useful..

这个报错主要是对称性设置,在`INCAR`中加入
```shell
ISYM = 0
SYMPREC = 1E-8
```
问题得到解决.

# 参考


