---
title: Latex 排版中的技巧记录(....ing)
tags: Study Latex
layout: article
license: true
toc: true
key: a20210130
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
这里是想通过这个博客来整理一下平时在使用Latex进行排版的时候会用到的一些特殊的方法,所以这个博客我用在自己使用Latex过程中,慢慢进行整理.
{:.info}
<!--more-->
# 前言
平时在使用Latex排版的时候总会遇到一些特殊的排版需求,锁着就想在这里以一个长期笔记的形式来将机子遇到的这些排版方法都整理到一起,也方便自己之后再时用的时候可以反过来看看.
# 公式图片并排
```latex
\begin{frame}
\frametitle{PRL,124,227001}
\begin{equation}
\begin{aligned}
H(\mathbf{k})&=2\lambda_x\sin k_x\sigma_xs_z\tau_z+2\lambda_y\sin k_y\sigma_y\tau_z\\
&+\left[(m_0-2t_x\cos k_x-2t_y\cos k_y)\sigma_z-\mu\right]\tau_z+\Delta_0\tau_x+\mathbf{h}\cdot\mathbf{s}
\end{aligned}
\end{equation}
\begin{columns}
	\column{0.5\textwidth}
	\begin{equation}
		H_{\mathrm{edge,j}}=-i\lambda_js_z\tau_z\partial_{l_j}+\Delta_0\tau_x+h_js_x
		\end{equation}
%After a unitary transformation $U=1\oplus(-is_y)$
		\begin{equation}
		H^{'}_{\mathrm{edge,j}}=-i\lambda_js_z\partial_{l_j}+\Delta_0s_x\tau_z+h_js_x
		\end{equation}
	\bigskip
	\column{0.4\textwidth}
	\begin{figure}[h]
	\centering
	\includegraphics[scale=0.15]{pic/c9.png}
\end{figure}
\end{columns}
	\pause
\begin{block}{}
{\color{blue!50!green}在这个模型中,利用Zeeman长辅助$s$-波超导体,在$\tau$空间中,在相邻的边界上也是形成了相反的"质量",将$\mathbf{h}\cdot\mathbf{s}$总体当作有效质量,其在BZ中的分布如上图.}
\end{block}
\end{frame}
```
**实现公式与图片放置到同一行的位置上. **

![png](/assets/images/latex/s1.png)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)