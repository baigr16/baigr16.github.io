---
title: Latex参考文献引用设置及补充材料公式编号修改
tags: Latex Study
layout: article
license: true
toc: true
pageview: true
key: a20200707
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
在科研写作中肯定要用到[Latex](https://www.latex-project.org/)的，正好最近遇到文献引用的问题，就想把自己使用的文献引用方式总结一些，并且探索了一下将**自己想解释的内容直接放到参考文献中并引用**，这种引用情况在阅读文献中经常看到，比如作者把某些内容放到了补充材料中，那么有时候在参考文献中就会看到写一段话，然后这段话解释了要去补充材料中哪里去寻找。
{:.success}
<!--more-->
# 文献引用
我通常使用的文献引用方式，是将参考文献单独写在一个文件中**(file.bib)**，里面的每条参考文献都会有唯一的引用号，这个引用号可以自行设定，正文中想要引用某一篇文献的时候，只需要**\cite{文献号}**即可，个人觉得这个方式很简单，而且**file.bib**中的内容可以直接从网站上面导出，下面以[PHYSICAL REVIEW B](https://journals.aps.org/prb/)为例，展示如何导出latex可用的文献引用。如下图所示，找到对应文献的导出连接，这里默认导出的是latex参考引用，点进去后也可以选择[EndNote](https://endnote.com/)的格式。

![png](/assets/images/latex/p1.png){:width="330px",:height="495px"}![png](/assets/images/latex/p2.png){:width="330px",:height="495px"}

以这种方式组织的bib文件，在进行参考文献索引时每个文献的标示号就是@article后的第一个参数，如图中所示这篇文章的标识号为*PhysRevB.102.020501*，所以在正文中想要引用这篇文献只需**\cite{PhysRevB.102.020501}**即可，但这里有个前提条件，想要使用bib文件来作为正文的引用，需要先在正文中说明，即在最后加入**\bibliography{ref}**，这里的ref就是你收集上面导出的参考文献信息的文件名。下面展示一个很简单的示例。
```latex
\documentclass[reprint,amsmath,amssymb,aps]{revtex4-2}
%========================package=====================
\usepackage{graphicx}% Include figure files
\usepackage{dcolumn}% Align table columns on decimal point
\usepackage{bm}% bold math
\usepackage[colorlinks=true,linkcolor=blue,anchorcolor=blue, citecolor=blue,urlcolor=blue]{hyperref}% add hypertext capabilities
\usepackage[mathlines]{lineno}% Enable numbering of text and display math
%-----------------------------------正文---------------------------------------------
\begin{document}
\section{sec1}
Hello\cite{PhysRevB.102.020501}
%-------------------------------
\bibliography{ref}  % 参考文献
\end{document}
```
这里就直接使用上图所示导出的一篇参考文献，上面就是这种方式参考文献是如何引用的，这里要保证**ref.bib**和你的**latex主文件**是在同一个文件夹下面，一般编译两次后就可以看到参考文献的引用了。(上面的latex代码不一定可以正确运行，我只是示例了一下这种方式中的参考文献是怎么引用，怎么组织的)

# Arxiv中文章的引用
如果你搞科研，肯定会遇到要引用Arxiv上面的文章，虽然它也提供了类似PRB一样的文献引用导出，它同样也有支持latex引用的文献格式(上面那种引用格式)，但是最后生成的参考文献的格式，和PRB上导出的那种，总会有一些差别，不会是你想要显示的方式。这里完全以PRB的参考文献格式来把Arxiv上的文章的参考引用，也和PRB一样。首先注意到的是，PRB导出的参考引用都是**@article**，虽然Arxiv也能导出类似的引用，但真的会遇到奇怪的问题，所以我索性放弃了Arxiv上的导出，直接自己写，但是这并不复杂，而且内容很少。这里使用**@misc**来整理Arxiv上文章的参考引用，而所需要整理的只有author、title、year、Eprint，在Eprint中的就是文章在Arxiv中对应的一个标示号，这是网站分配给每一个预印本文章的。示例如下
```latex
@misc{re127,
Author = {Maryam Khezerlou and Hadi Goudarzi},
Title = {Mass-like gap creation by mixed singlet and triplet state in superconducting topological insulator},
Year = {2019},
Eprint = {arXiv:1903.01144},
}
```
通过上面这种方式就可以将文献的参考引用格式整理成和PRB导出的格式类似的形式，这里要说明一下，按照这种形式引用的Arxiv文章，年份会在引用文献title的后面，而标准的PRB参考文献的年份是在期刊名后面的，但是我看到PRL上文章的参考引用中，Arxiv上的文章并没有给出年份，所以这里可以直接把year删除了，也是可以的。
# 任意内容的引用
在开头的时候提到过作者有时候会写一段话来作为参考引用，这中引用方式其实也是通过**@misc**来实现，只不过这时候需要输入的内容就只有你想解释的文字而已。
```latex
@misc{re3,
  title={This is....},
}
```
在ref.bib文件中加入这种形式的参考后，你就可以在正文想要引用的地方通过**\cite{re3}**来直接引用这段话。到这里我自己关于参考文献引用的探索就结束了，如果之后有新的内容，还会继续补充。

# 自定义公式编号

最近在整理文章的补充材料，正好遇到了要修改公式编号的问题。通常正文中的公式编号都是按照公式的顺序号来的，第几个公式就是第几号(**我这里是以revtex4为模板说的，当然你也可以让公式编号按照章节走，这是写毕业论文的格式**)，但是如果你加**\appendix**，那么附录中每一节(section)的公式编号会有个大写字母来表示(A1,B1，C1等等)，但是到了补充材料中，通常看到公式都是以S为开头，之后就是公式的顺序编号了，所以这里可以人为的修改这个公式编号，来满足需求，Latex真的是太方便了。

```latex
\renewcommand{\theequation}{S\arabic{equation}} % This line ads "S" in front of your equation numbering.
```

在导言区加入这句话之后，你之后的所有公式都以S开头的。这里要注意，这个和附录不要同时使用，不然的话附录中的每个section中的公式编号会是独立的，每个section的公式编号会以S开头，但是都会重新以1为起始重新编号。上面的设置参考了[这里]( https://tex.stackexchange.com/questions/164640/customize-equation-numbering-for-equation-environment ).

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)