---
title: 哈密顿量降维的方法
tags:  Study
layout: article
license: true
toc: true
key: a20220301
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
这里整理一下如何将一个$4\times 4$或者更大的哈密顿量分解成$2\times 2$的维度更低的哈密顿量来对哈密顿量的一些性质进行研究,这个方法在文献中是比较常见的一种处理方式.
{:.info}
<!--more-->
# 问题说明
在研究有些问题的时候, 碰到的哈密顿量假如是$4\times 4$的形式, 比如下面的哈密顿量

$$
\begin{equation}
H=v_fk\sigma_3\tau_3+h\sigma_1+\Delta\tau_1\label{q1}
\end{equation}
$$

简单的描述一下,这其实就是实现拓扑超导的一个模型,$h$表示磁场,$\Delta$表示电子配对,这个哈密顿量就是$4\times 4$的,而且从对易关系上来看,这两项之间满足对易关系,也就说对于边界态$v_fk\sigma_3\tau_3$而言,这里的两个质量项$h\sigma_1$和$\Delta\tau_1$之间是竞争关系,即满足

$$
[\sigma_1\otimes\tau_0,\sigma_0\otimes \tau_1]=0
$$

这时候总的质量项就会是

$$
h\pm\Delta
$$

这里其实都是根据质量项之间的对易关系来分析的,如果熟悉对易关系和反对易关系的话,上面的结论上比较容易得到的,但是这里如果可以从一个更加数学的角度来得到这个结论就会更加的清晰.

首先来回顾一下量子力学,如果找到一个量$O$和哈密顿量之间满足对易关系

$$
[H,O]=0
$$

这个时候两个算符具有共同本征态,就可以利用$O$的本征值来标记哈密顿量$H$的本征态,比如算符$O$是反演操作算符,那么就可以对$H$的每个本征态标记其宇称$\pm\xi$,这里的宇称也就是这个态对应的反演算符的本征值.

有了上面的准备,这个时候就可以对哈密顿量(\ref{q1})也来寻找到一个和其对易的项,来对这个哈密顿量的本征态进行分类.首先再来看看哈密顿量

$$
\begin{equation}
H=v_fk\sigma_3\tau_3+h\sigma_1+\Delta\tau_1
\end{equation}
$$

说白了就是用一堆Pauli矩阵将哈密顿量表示了出来,那么在这里就可以找一个Pauli矩阵和哈密顿量满足对弈关系

$$
[\Gamma_i,H]=0
$$

因为上面的哈密顿量(\ref{q1})形式简单,可以发现$\Gamma_i=\sigma_1\tau_1$就满足这样的关系,而且$\Gamma_i$的本征值为$\pm 1$,所以就可以分别得到这个$\Gamma_i$矩阵的本征态了,那么接下来就将哈密顿量分别投影到本征值对应的本征空间中,假设$\Gamma_i$的本征态标记为
$\rvert\psi^\pm\rangle$,可以得到对于本征值为$\pm 1$的本征态,对应的哈密顿量为

$$
\begin{equation}
\begin{aligned}
H^+&=\langle\psi^+\rvert H\rvert\psi^+\rangle=v_fk\tilde{\sigma}_3+(h+\Delta)\tilde{\sigma}_2\\
H^-&=\langle\psi^-\rvert H\rvert\psi^-\rangle=v_fk\tilde{\sigma}_3+(h-\Delta)\tilde{\sigma}_2
\end{aligned}
\end{equation}
$$
{:.success}

到这里就可以看到在这种形式下就将原本$4\times 4$的哈密顿量分解成了两个$2\times 2$的哈密顿量,而且质量项分别为

$$
(h+\Delta)\tilde{\sigma}_2\qquad (h-\Delta)\tilde{\sigma}_2
$$

这里就从数学上将前面直接利用对易关系来分析质量项之间的关系明确的表示出来.

# 代码
这里截图给一下利用Mathematica计算的结果

![png](/assets/images/Majorana/f2.png)

![png](/assets/images/Majorana/f1.png)

# 参考
 - [Classification of topological quantum matter with symmetries](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.88.035005)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)