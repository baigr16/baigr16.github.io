---
title: 使用Latex进行文档对比并自动标记修改的内容
tags:  Latex
layout: article
license: true
toc: true
key: a20240331
pageview: true
# cover: /assets/images/Julia/julia-logo.png
header:
  theme: dark
  background: 'linear-gradient(to right, rgb(101, 78, 163), rgb(234, 175, 200))'
article_header:
  type: overlay
  theme: dark
  background_color: false
  background_image: 
    gradient: 'linear-gradient(to right, #d3cce3, #e9e4f0)'
    image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
这里整理一下如何使用Latex自动整理出文档改动内容并进行标记，这个功能在我们给审稿人回复意见的时候非常有用。
{:.info}
<!--more-->
# 前言
通常在文章审稿的过程中，都需要根据审稿人的要求和问题对文章内容进行修改，常用的方法就是将修改的部分用颜色标记出来，但是这样就将原来的内容删除掉了，审稿人再次看到的时候只能看到我们修改文章的结果，但具体修改了什么，修改到了什么程度就不是很显而易见了。当然，如果有精力也可以用下划线以及删除线等标记自己手动将原本的内容与修改之后的内容分别进行标记，但这种手动的方式操作起来还是挺累的，毕竟`Latex`本身虽然排本能力很强，但是在写的时候可读性就不是很好了，内容多了就会眼花缭乱。

# 解决方法
实际上在安装了`TexLive`之后，本身就会自带一个文档比对的功能，与`Linux`系统里面的`diff`的功能是相似的，使用
```shell
latexdiff file-1.tex file-2.tex > diff.tex
```
这个命令之后，就可以自动生成一个对比之后的文件`diff.tex`，当然，这里的名字都是自己起的。在产生的`diff.tex`文件中就会将两个文件的比对结果存储，删除以及修改的内容都会用各种标记方式给出。再对`diff.tex`文件编译之后就可以得到一个内容修改对比的结果了，示例如下

- 原版(file-1.tex)

![](https://files.mdnice.com/user/42079/f69d4bd5-3486-469f-8b3d-87a04371db49.png)

- 修改版(file-1.tex)

![](https://files.mdnice.com/user/42079/f79e6f35-40a1-4f91-becd-ae76fda7b221.png)

- 对比版(diff.tex)
  
![](https://files.mdnice.com/user/42079/b69df5af-5ecb-47a2-bf9f-afbd760fdba3.png)

这样就可以让审稿人一目了然的看清楚我们对正文的修改以及修改程度，对文章审稿意见的仔细回复也能体现出对审稿人的尊重
{:.info}

上面给出的只是`latexdiff`的默认参数选择，更加丰富的选项可以移步官网查看，但是以我现在的需要，看起来默认的选项就已经足够了。

# 提示
前面只是给出了一个很简单的示例，实际上文档中包含的并不仅仅是文字，还会有公式和图片，这些都是可以识别修改前后的差别。

但有一点很重要，在修改之后的*file-2.tex*中，尽量要保证它与*file-1.tex*的结构是一致的。比如说你有一张图片本来在Latex中在*Section I*里面的，但是在修改过程中将这个图片的插入移动到了*Section II*中，此时再利用上面的方法进行文档差异识别的时候，就会将*file-2.tex*中的这张图片插入与*file-1.tex*中处于相同文本位置的文字进行对比，给出的结果自然是有问题的。因此在修改的时候，若非文章进行的改动非常大，尽量不要去改变原本的架构。
{:.warning}



# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg){:.border.rounded}{:width="300px" height="300px"}
<div class="card">
  <div class="card__content">
    <div class="card__header">
      <h4>Email</h4>
    </div>
    <p>yxli406@gmail.com</p>
  </div>
</div>