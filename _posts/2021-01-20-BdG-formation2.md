---
title: 超导态基矢选择对构建BdG哈密顿量的影响
tags:  Superconductor
layout: article
license: true
toc: true
key: a20210120b
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
在[正常态到BdG哈密顿量的构建](https://yxli8023.github.io/2021/01/20/BdG-formation.html)这篇博客中我介绍了如何将正常态的哈密顿量在加入超导配对之后改写为BdG的形式,这里来讨论另外一种常见的情形,那就是基矢选择对于这个构建过程有什么影响.
{:.info}
<!--more-->
# 正常态模型
这里选择正常态哈密顿量为

$$h^{\mathrm{TI}}(\mathbf{k})=\left[2t(\cos k_x-\cos k_y)+4t_1\cos k_x \cos k_y\right]\sigma_z+2\lambda(\sin k_xs_y-\sin k_ys_x)\sigma_x\label{nor}$$

take $m(\mathbf{k})=2t(\cos k_x-\cos k_y)+4t_1\cos k_x\cos k_y$.

这个模型同样是一个量子自旋效应的模型,不过和[BHZ模型](https://topocondmat.org/w6_3dti/bhz.html)不同的是,它的边界态对于不同情况下的开边界条件,出现的位置是不同的,这篇博客中[正常态到BdG哈密顿量的构建](https://yxli8023.github.io/2021/01/20/BdG-formation.html)的模型,无论哪个方向取为开边界,边界态出现的位置都是相同的,具体的内容可以查看**参考中的(3),(4)**两篇文章.

这里就用和[正常态到BdG哈密顿量的构建](https://yxli8023.github.io/2021/01/20/BdG-formation.html)这篇博客中相同的方法来推导$(\ref{nor})$在加入超导配对之后的BdG形式

$$\begin{equation}
\begin{aligned}
\sum_\mathbf{k}m(\mathbf{k})c^\dagger_{a\uparrow\mathbf{k}}c_{a\uparrow\mathbf{k}}&=\frac{1}{2}\sum_\mathbf{k}m(\mathbf{k})\left[c^\dagger_{a\uparrow\mathbf{k}}c_{a\uparrow\mathbf{k}}-c_{a\uparrow\mathbf{k}}c^\dagger_{a\uparrow\mathbf{k}}\right]\\
&=\frac{1}{2}\sum_\mathbf{k}\left[m(\mathbf{k})c^\dagger_{a\uparrow\mathbf{k}}c_{a\uparrow\mathbf{k}}-m(\mathbf{-k})c_{a\uparrow\mathbf{-k}}c^\dagger_{a\uparrow\mathbf{-k}}\right]\\
&=\frac{1}{2}\sum_\mathbf{k}m(\mathbf{k})\left[c^\dagger_{a\uparrow\mathbf{k}}c_{a\uparrow\mathbf{k}}-c_{a\uparrow\mathbf{-k}}c^\dagger_{a\uparrow\mathbf{-k}}\right]
\end{aligned}
\end{equation}\label{m21}$$

$$\begin{equation}
\begin{aligned}
\sum_\mathbf{k}-i\sin k_x c^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}&=\frac{1}{2}\sum_\mathbf{k}-i\sin k_x\left[c^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}-c_{b\downarrow\mathbf{k}}c^\dagger_{a\uparrow\mathbf{k}}\right]\\
&=\frac{1}{2}\sum_\mathbf{k}\left[-i\sin k_xc^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}+i\sin (-k_x)c_{b\downarrow\mathbf{-k}}c^\dagger_{a\uparrow\mathbf{-k}}\right]\\
&=\frac{1}{2}\sum_\mathbf{k}-i\sin k_x\left[c^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}-c_{b\downarrow\mathbf{-k}}c^\dagger_{a\uparrow\mathbf{-k}}\right]
\end{aligned}
\end{equation}\label{m22}$$

$$\begin{equation}
\begin{aligned}
-\sum_\mathbf{k}\sin k_yc^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}&=-\frac{1}{2}\sum_\mathbf{k}\sin k_y\left[c^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}-c_{b\downarrow\mathbf{k}}c^\dagger_{a\uparrow\mathbf{k}}\right]\\
&=-\frac{1}{2}\sum_\mathbf{k}\left[\sin k_yc^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}-\sin(-k_y)c_{b\downarrow\mathbf{-k}}c^\dagger_{a\uparrow\mathbf{-k}}\right]\\
&=-\frac{1}{2}\sum_\mathbf{k}\sin k_y\left[c^\dagger_{a\uparrow\mathbf{k}}c_{b\downarrow\mathbf{k}}+c_{b\downarrow\mathbf{-k}}c^\dagger_{a\uparrow\mathbf{-k}}\right]
\end{aligned}
\end{equation}\label{m23}$$

$$\begin{equation}
\begin{aligned}
\sum_\mathbf{k\alpha}\Delta(\mathbf{k})c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}&=\frac{1}{2}\sum_{\mathbf{k}\alpha}\Delta(\mathbf{k})\left[c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}-c^\dagger_{\alpha\downarrow\mathbf{-k}}c^\dagger_{\alpha\uparrow\mathbf{k}}\right]\\
&=\frac{1}{2}\sum_{\mathbf{k}\alpha}\left[\Delta(\mathbf{k})c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}-\Delta(\mathbf{-k})c^\dagger_{\alpha\downarrow\mathbf{k}}c^\dagger_{\alpha\uparrow\mathbf{-k}}\right]\\
&=\frac{1}{2}\sum_{\mathbf{k}\alpha}\Delta(\mathbf{k})\left[c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}-c^\dagger_{\alpha\downarrow\mathbf{k}}c^\dagger_{\alpha\uparrow\mathbf{-k}}\right]
\end{aligned}
\end{equation}\label{m24}$$

$\Gamma_1=\tau_z\otimes s_0\otimes\sigma_z\quad\Gamma_2=\tau_z\otimes s_y\otimes\sigma_x\quad\Gamma_3=\tau_z\otimes s_x\otimes\sigma_x\quad\Gamma_4=\tau_z\otimes s_0\otimes \sigma_0$

$$\begin{equation}
\begin{aligned}
(\ref{m21})=\Psi^\dagger\Gamma_1\Psi\quad(\ref{m22})=\Psi^\dagger\Gamma_2\Psi\quad(\ref{m23})=\Psi^\dagger\Gamma_3\Psi\quad(\ref{m24})=\Psi^\dagger\Gamma_4\Psi
\end{aligned}
\end{equation}$$

$$\begin{equation}
\begin{aligned}
H^{\mathrm{BdG}}(\mathbf{k})=&(h^{\mathrm{TI}}(\mathbf{k})-\mu)\tau_z+\Delta(\mathbf{k})\tau_x\\
h^{\mathrm{TI}}(\mathbf{k})=&\left[2t(\cos k_x-\cos k_y)+4t_1\cos k_x \cos k_y\right]\sigma_z+2\lambda(\sin k_xs_y-\sin k_ys_x)\sigma_x\\
\Delta(\mathbf{k})=&\Delta_0+2\Delta_1(\cos k_x+\cos k_y)
\end{aligned}
\end{equation}$$

Basis for $H^{\mathrm{BdG}}(\mathbf{k})$ is  $\Psi^\dagger=(c^\dagger_{a\uparrow\mathbf{k}},c^\dagger_{b\uparrow\mathbf{k}},c^\dagger_{a\downarrow\mathbf{k}},c^\dagger_{b\downarrow\mathbf{k}},c_{a\downarrow\mathbf{-k}},c_{b\downarrow\mathbf{-k}},-c_{a\uparrow\mathbf{-k}},-c_{b\uparrow\mathbf{-k}})=(C_\mathbf{k}^\dagger,-is_y\sigma_0C_\mathbf{-k})$, with $C^\dagger_\mathbf{k}=(c^\dagger_{a\uparrow\mathbf{k}},c^\dagger_{b\uparrow\mathbf{k}},c^\dagger_{a\downarrow\mathbf{k}},c^\dagger_{b\downarrow\mathbf{k}})$

# 过程分析
这里的推导过程和[正常态到BdG哈密顿量的构建](https://yxli8023.github.io/2021/01/20/BdG-formation.html)这篇博客完全相同,不过要注意的是这里的基矢选择是和前面不同的

Basis for $H^{\mathrm{BdG}}(\mathbf{k})$ is  $\Psi^\dagger=(c^\dagger_{a\uparrow\mathbf{k}},c^\dagger_{b\uparrow\mathbf{k}},c^\dagger_{a\downarrow\mathbf{k}},c^\dagger_{b\downarrow\mathbf{k}},c_{a\downarrow\mathbf{-k}},c_{b\downarrow\mathbf{-k}},-c_{a\uparrow\mathbf{-k}},-c_{b\uparrow\mathbf{-k}})=(C_\mathbf{k}^\dagger,-is_y\sigma_0C_\mathbf{-k})$, with $C^\dagger_\mathbf{k}=(c^\dagger_{a\uparrow\mathbf{k}},c^\dagger_{b\uparrow\mathbf{k}},c^\dagger_{a\downarrow\mathbf{k}},c^\dagger_{b\downarrow\mathbf{k}})$
{:.warning}
正式因为这个基矢选择的不同,所以最终的结果就是

$$H^{\mathrm{BdG}}(\mathbf{k})=(h^{\mathrm{TI}}(\mathbf{k})-\mu)\tau_z+\Delta(\mathbf{k})\tau_x$$

在这个基矢下,正常态和超导态看起来会比较直接,$\tau$代表的粒子空穴自由度,这时候正常态相当于处在$\tau_z$这个表示下,而$\Delta(\mathbf{k})$则是处在$\tau_x$这个表示下,相当于超导这一项就是将正常态的粒子部分和空穴部分耦合起来,因为将它们完整写出来,$\Delta(\mathbf{k})$正好就是在非对角线上.

$$\begin{equation}\left[\begin{array}{cc}h^{\mathrm{TI}}(\mathbf{k})-\mu&\Delta(\mathbf{k})\\\Delta^\dagger(\mathbf{k})&-(h^{\mathrm{TI}}(\mathbf{k})-\mu)\end{array}\right]\end{equation}$$

在这个表示下,超导配对和正常态看起来就是很直观的,但是在写矩阵的时候,不仅需要考虑哈密顿量中的矩阵,还需要考虑基矢中的符号,这个问题在以前的博客中有提及到,可以参考这里[Hamiltonian构建时的基矢选择](https://yxli8023.github.io/2020/07/03/Basis-Chose.html)这篇博客,我在这里曾指出过这个问题.

# 参考
1.[Bogoliubov变换与Majorana表示](https://zhuanlan.zhihu.com/p/59445571)

2.[Bogoliubov-de Gennes Method and Its Applications](https://link.springer.com/book/10.1007/978-3-319-31314-6)

3.[High-Temperature Majorana Corner States](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.186801)

4.[Majorana Corner Modes in a High-Temperature Platform](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.096803)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)