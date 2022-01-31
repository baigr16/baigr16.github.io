---
title: 正常态到BdG哈密顿量的构建
tags:  Superconductor
layout: article
license: true
toc: true
key: a20210120a
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
在平均场的方法下研究叉掉问题,通常需要将哈密顿量写成[BdG哈密顿量](https://link.springer.com/book/10.1007/978-3-319-31314-6)的形式,其实也就是用时考虑电子和空穴算符,这样就会使得原来的哈密顿量矩阵扩大一倍.这里就详细的从一个模型出发,推导一下如何将一个正常态的哈密顿量,在加入超导配对后改写成BdG形式的哈密顿量.
{:.info}
<!--more-->
# BHZ模型
这里利用拓扑绝缘体的[BHZ模型](https://topocondmat.org/w6_3dti/bhz.html)来演示如何将它在加入超导之后变成一个BdG的形式,首先BHZ模型的哈密顿量为

$$H_{\mathrm{TI}}=(m_0-t_x\cos k_x-t_y\cos k_x)\sigma_z+\lambda_x\sin k_x\sigma_x s_z+\lambda_y\sin k_y\sigma_y-\mu\label{ti}$$

with basis $\Psi=(c_{a\uparrow\mathbf{k}},c_{b\uparrow\mathbf{k}},c_{a\downarrow\mathbf{k}},c_{b\downarrow\mathbf{k}})^T$, $M(\mathbf{k})=m_0-t_x\cos k_x-t_y\cos k_y$

将哈密顿量以算符形式写出来即为

$$\hat{H}^{TI}=\sum_\mathbf{k}\Psi^\dagger H^{\mathrm{TI}}(\mathbf{k})\Psi$$

当在模型中加入超导配对之后,模型变为

$$H^{\mathrm{BdG}}(\mathbf{k})=(m_0-t_x\cos k_x-t_y\cos k_y)\sigma_z\tau_z+\lambda_x\sin k_x\sigma_xs_z+\lambda_y\sin k_y\sigma_y\tau_z+\Delta(\mathbf{k})s_y\tau_y-\mu\tau_z\label{bdg}$$

with basis $\Psi=(c_{a\uparrow\mathbf{k}},c_{b\uparrow\mathbf{k}},c_{a\downarrow\mathbf{k}},c_{b\downarrow\mathbf{k}},c^\dagger_{a\uparrow\mathbf{-k}},c^\dagger_{b\uparrow\mathbf{-k}},c^\dagger_{a\downarrow\mathbf{-k}},c^\dagger_{b\downarrow\mathbf{-k}})^T=(C_\mathbf{k},C^\dagger_\mathbf{-k})$

将哈密顿量以算符形式写出来即为

$$\hat{H}^{BdG}=\frac{1}{2}\sum_\mathbf{k}\Psi^\dagger H^{\mathrm{BdG}}(\mathbf{k})\Psi$$

下面就来推导每一项到底是如何从公式$(\ref{ti})$在加入超导之后变为$(\ref{bdg})$的.

$$\begin{eqnarray}
\begin{aligned}
\sum_\mathbf{k}M(\mathbf{k})c^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}&=\frac{1}{2}\sum_\mathbf{k}\left[M(\mathbf{k})c^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}-M(\mathbf{k})c_{b\sigma\mathbf{k}}c^\dagger_{a\sigma\mathbf{k}}\right]\\
&=\frac{1}{2}\sum_\mathbf{k}\left[M(\mathbf{k})c^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}-M(\mathbf{-k})c_{b\sigma\mathbf{-k}}c^\dagger_{a\sigma\mathbf{-k}}\right]
\end{aligned}\label{q1}
\end{eqnarray}$$

$$\begin{equation}
\begin{aligned}
\sum_\mathbf{k}\sin k_xc^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}&=\frac{1}{2}\sum_\mathbf{k}\left[\sin k_xc^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}-\sin k_xc_{b\sigma\mathbf{k}}c^\dagger_{a\sigma\mathbf{k}}\right]\\
&=\frac{1}{2}\sum_\mathbf{k}\left[\sin k_xc_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}-\sin (-k_x)c_{b\sigma\mathbf{-k}}c^\dagger_{a\sigma\mathbf{-k}}\right]\\
&=\frac{1}{2}\sum_\mathbf{k}\left[\sin k_xc_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}+\sin (k_x)c_{b\sigma\mathbf{-k}}c^\dagger_{a\sigma\mathbf{-k}}\right]\label{m12}
\end{aligned}
\end{equation}\label{q2}$$

$$\begin{equation}
\begin{aligned}
\sum_{\mathbf{k}\sigma}-i\sin k_yc^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}&=\frac{1}{2}\sum_{\mathbf{k}\sigma}\left[-i\sin k_yc^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}+i\sin k_yc_{b\sigma\mathbf{k}}c^\dagger_{a\sigma\mathbf{k}}\right]\\
&=\frac{1}{2}\sum_{\mathbf{k}\sigma}\left[-i\sin k_xc^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}+i\sin(-k_y)c_{b\sigma\mathbf{-k}}c^\dagger_{a\sigma\mathbf{-k}}\right]\\
&=\frac{1}{2}\sum_{\mathbf{k}\sigma}\left[-i\sin k_xc^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}-i\sin(k_y)c_{b\sigma\mathbf{-k}}c^\dagger_{a\sigma\mathbf{-k}}\right]\label{m13}
\end{aligned}
\end{equation}\label{q3}$$

$$\begin{equation}
\begin{aligned}
\sum_\mathbf{k\alpha}\Delta(\mathbf{k})c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}&=\frac{1}{2}\sum_\mathbf{k\alpha}\left[\Delta(\mathbf{k})c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}-\Delta(\mathbf{k})c^\dagger_{\alpha\downarrow\mathbf{-k}}c^\dagger_{\alpha\uparrow\mathbf{k}}\right]\\
&=\frac{1}{2}\sum_{\mathbf{k}\sigma}\left[\Delta(\mathbf{k})c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}-\Delta(\mathbf{-k})c^\dagger_{\alpha\downarrow\mathbf{k}}c^\dagger_{\alpha\uparrow\mathbf{-k}}\right]\\
&=\frac{1}{2}\sum_{\mathbf{k}\sigma}\left[\Delta(\mathbf{k})c^\dagger_{\alpha\uparrow\mathbf{k}}c^\dagger_{\alpha\downarrow\mathbf{-k}}-\Delta(\mathbf{k})c^\dagger_{\alpha\downarrow\mathbf{k}}c^\dagger_{\alpha\uparrow\mathbf{-k}}\right]\label{m14}
\end{aligned}
\end{equation}\label{q4}$$

$\Gamma_1=\tau_z\otimes s_0\otimes\sigma_z\quad\Gamma_2=\tau_0\otimes s_z\sigma_x\quad\Gamma_3=\tau_z\otimes s_0\otimes\sigma_y\quad\Gamma_4=\tau_y\otimes s_y\otimes\sigma_0$

$$(\ref{q1})=\Psi^\dagger\Gamma_1\Psi\quad(\ref{q2})=\Psi^\dagger\Gamma_2\Psi\quad(\ref{q3})=\Psi^\dagger\Gamma_3\Psi\quad(\ref{q4})=\Psi^\dagger\Gamma_4\Psi$$

# 过程分析
这里从正常态加入超导配对之后最主要的变化就是将正常态的项,拆分为粒子部分$M(\mathbf{k})c^\dagger_{a\sigma\mathbf{k}}c_{b\sigma\mathbf{k}}$和空穴部分$M(\mathbf{k})c_{b\sigma\mathbf{k}}c^\dagger_{a\sigma\mathbf{k}}=M(\mathbf{-k})c_{b\sigma\mathbf{-k}}c^\dagger_{a\sigma\mathbf{-k}}$,也就是说将粒子部分动量$\mathbf{k}$,单独分出来然后变成对应的空穴部分$\mathbf{-k}$,最后会出现一个$1/2$的系数,但是作为哈密顿量,乘以一个常数,或者加一个常数,对本征的结果是没有影响的,所以如果你不追求严格性,将多余的常数丢弃,将前面的$1/2$这个系数不屑,结果也是正确的,因为这个只会影响你取势能零点的位置,对绝对的能量则是没有影响的.

其实在公式$(\ref{q1})$中,然全是利用的费米子的反对易关系才可以将算符分解成那样的形式
$$\{c_k^\dagger,c_k\}=1\rightarrow c^\dagger_kc_k+c_kc_k^\dagger=1\rightarrow c^\dagger_kc_k=1-c_kc_k^\dagger=1-c_{-k}c_{-k}^\dagger$$
这里我直接将常数项*1*扔去,没有写出来,所以严格的形式应该为
$$\hat{H}^{BdG}=\frac{1}{2}\sum_\mathbf{k}\Psi^\dagger H^{\mathrm{BdG}}(\mathbf{k})\Psi+\mathrm{const}$$

以上就是如何将一个正常态哈密顿量在加入超导配对之后,写成BdG形式的全部过程.

# 参考
1.[Bogoliubov变换与Majorana表示](https://zhuanlan.zhihu.com/p/59445571)

2.[Bogoliubov-de Gennes Method and Its Applications](https://link.springer.com/book/10.1007/978-3-319-31314-6)

3.[Majorana Corner Modes in a High-Temperature Platform](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.096803)