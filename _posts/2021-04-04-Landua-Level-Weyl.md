---
title: Weyl半金属中朗道能级求解
tags: Topology 
layout: article
license: true
toc: true
key: a20210404
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
pageview: true
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
这里整理了在Weyl半金属中加入磁场之后,哈密顿量的朗道能级如何求解,只是利用了最常用的解方程的方法.
<!--more-->
自由电子气中加入磁场会产生朗道量子化能级,同样的对于其他体系也会有这样的情况发生,这里就整理一下最近给别人回答问题的时候,计算了Weyl半金属加入磁场之后的朗道能级求解.
{:.info}

这里想要求解的内容来自于参考文献1,我暂时不清楚文章中利用算符变换来求解得到朗道能级的方法,我是从解方程的角度直接求解的.
# 求解
首先加入磁场之后哈密顿量要进行变换

$$H_\chi({\bf k})=\chi v_F(\hbar{\bf k}+e{\bf A})/cdot \sigma$$

磁场${\bf B}=(B,0,0)$, 且${\bf B} = \nabla\times{\bf A}$, 这里将矢势${\bf A}=(0,0,By)$, 则哈密顿量为

$$H_\chi({\bf k})=\chi v_F\left[\hbar k_x\sigma_x+\hbar k_y\sigma_y+\hbar(k_z-\frac{eB}{\hbar}y)\sigma_z\right]\label{e1}$$

此时来求解(\ref{e1})对应的本征方程

$$H_\chi({\bf k})\psi=E\psi$$

这里我不去求解关于$H$的本征方程, 而是将$H$再作用一次, 求解$H^2\psi=E^2\psi$的本征方程.
{:.warning}

$$\begin{equation}\begin{aligned}H_\chi^2({\bf k})&=\chi^2v_F^2\left[\hbar^2k_x^2+\hbar^2k_y^2+\hbar^2(k_z-\frac{eB}{\hbar}y)^2\right]\\
H_\chi^2({\bf k})\psi&=\chi^2v_F^2\left[\hbar^2k_x^2+\hbar^2k_y^2+\hbar^2(k_z-\frac{eB}{\hbar}y)^2\right]\psi=E^2\psi\end{aligned}\end{equation}$$

$$\left[\chi^2v_F^2\hbar^2k_y^2+\chi^2\hbar^2v_F^2(\frac{eB}{\hbar})^2(\frac{\hbar}{eB}k_z-y)^2\right]=(E^2-\chi^2v_F^2\hbar^2k_x^2)\psi=\epsilon\psi$$

化简到这一步, 可以看到和自由电子气模型加入磁场的情况变得完全相同, 那么将这个形式进一步修正, 变成自由电子气模型

$$\left[-\frac{\hbar^2}{\hbar/\chi^2v_F^2}\frac{\partial^2}{\partial y^2}+\frac{\frac{\hbar}{2\chi^2v_F^2}}{2}\cdot(2\frac{\chi^2v_F^2eB}{\hbar})^2(\frac{\hbar}{eB}k_z-y)^2\right]=\epsilon\psi$$

这里与自由电子气模型相比较$m=\frac{\hbar}{2\chi^2v_F^2}$, $\omega_0=\frac{2\chi^2v_F^2eB}{\hbar}$, $\epsilon=(n+\frac{1}{2})\omega_0$

最终可以求解得到能量为

$$E=\pm\sqrt{(n+\frac{1}{2})\hbar\omega_0+\chi^2v_F^2\hbar^2k_x^2}$$

原文中的结果为

$$\epsilon^\chi_n(k_x)=\left \{ \begin{array}{c}-\chi\hbar v_Fk_x\qquad n=0\\\text{sgn}(n)\sqrt{2\rvert n\rvert(\hbar\omega_c)^2+(\hbar v_Fk_x)^2}\qquad n\neq 0\end{array}\right.$$

这里$\omega_c=v_f/\mathcal{l}_B$, $\mathcal{l}_B=\sqrt{\hbar/eB}$. 

我这里求解得到的结果与原文中其实是完全相同的, 只不过原文中是利用算符变换求解得到的, 而且原文中应该是同时丢弃了零点能, 这样我的求解结果就是原文是完全一致的.
{:.warning}

# 升降算符法求解
在这里为了避免符号繁琐, 采用自然单位制$\hbar=e=1$, 并将手性$\xi$略去, 为了结果之后只需要简单在结果前面乘上这一项就可以, 磁场矢势此时采用对称规范$\mathbf{A}=\frac{1}{2}(0,-Bz,+By)$.

$$H(\mathbf{k})=k_x\sigma_x+(k_y-\frac{1}{2}Bz)\sigma_y+(k_z+\frac{1}{2}By)\sigma_z$$

重新定义新的算符$k_y^{'}=k_y-\frac{1}{2}Bz,k_z^{'}=k_z+\frac{1}{2}By$, 两者之间满足对易关系

$$\left[k_y^{'},k_z^{'}\right]=iB$$

接下来对哈密顿量进行一个幺正变换, $\sigma_x\rightarrow-\sigma_z,\sigma_y\rightarrow\sigma_x,\sigma_z\rightarrow-\sigma_y$, 利用到的幺正变换矩阵为

$$\begin{equation}U=\frac{1}{\sqrt{2}}\left[\begin{array}{cc}
i&-i\\
1&1
\end{array}\right]\end{equation}$$

哈密顿量改写为

$$H(\mathbf{k})=-k_x\sigma_z+k_y^{'}\sigma_x-k_z^{'}\sigma_y$$

定义升降算符

$$a=\frac{k_y^{'}-ik_z^{'}}{\sqrt{2B}}\qquad a^\dagger=\frac{k_y^{'}+ik_z^{'}}{\sqrt{2B}}$$

哈密顿量的矩阵形式为

$$H(\mathbf{k})=\sqrt{2B}\left(\begin{array}{cc}-\frac{k_x}{\sqrt{2B}}&a\\
a^\dagger&\frac{k_x}{\sqrt{2B}}\end{array}\right)$$

求解本征方程

$$H(\mathbf{k})\Psi=E\Psi\qquad \Psi=(\theta_1\rvert n-1\rangle,\theta_2\rvert n\rangle)^T$$

结合升降算符性质

$$a\rvert n\rangle=\sqrt{n}\rvert n-1\rangle\qquad a^\dagger\rvert n\rangle=\sqrt{n+1}\rvert n+1\rangle$$

可以求解得到本征值为

$$E=\pm\sqrt{k_x^2+2B n}\qquad (n\ge 1)$$

这里$n\ge 1$是因为波函数形式为$\Psi=(\theta_1\rvert n-1\rangle,\theta_2\rvert n\rangle)^T$, 所以必须满足这个条件. 若想要求解$n=0$时的能级, 则$\rvert n-1\rangle\rightarrow_{n=0}0$, 则可以求解得到

$$E=-k_x$$

最后求解得到的结果为

$$E(k_x)=\left\{\begin{array}{c}-k_x\qquad n=0\\
\sqrt{k_x^2+2B n}\qquad (n\ge 1)\end{array}\right.$$

# 矩阵方法
其实刚才已经将哈密顿量变换成矩阵算符形式

$$H(\mathbf{k})=\sqrt{2B}\left(\begin{array}{cc}-\frac{k_x}{\sqrt{2B}}&a\\
a^\dagger&\frac{k_x}{\sqrt{2B}}\end{array}\right)$$

把波函数写作

\Psi=(\theta_1\rvert n-1\rangle,\theta_2\rvert n\rangle)^T

我们可以发现升降算符只出现在非对角线上, 也就是$H(\mathbf{k})$作用到$\Psi$上的时候, $a^\dagger$作用到$\rvert n-1\rangle$, 而 $a$作用到了$\rvert n\rangle$上, 再结合升降算符的性质, 可以直接把矩阵写成

$$H(\mathbf{k})=\sqrt{2B}\left(\begin{array}{cc}-\frac{k_x}{\sqrt{2B}}&\sqrt{n}\\
\sqrt{n}&\frac{k_x}{\sqrt{2B}}\end{array}\right)$$

这个时候哈密顿量的本征值可以通过直接对角化这个$2\times 2$的矩阵得到, 结果和上面是完全相同的.

# 参考

1. [Quantum Oscillations of the Positive Longitudinal Magnetoconductivity: A Fingerprint for Identifying Weyl Semimetals](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.036601)

2. Topological insulator and topological superconductor

3.[Dirac Quantum Wells at Domain Walls in Antiferromagnetic Topological Insulators](https://arxiv.org/pdf/2104.00690.pdf)
