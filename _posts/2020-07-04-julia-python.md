---
title: Julia调用Python画图
tags: Study Julia Code 
layout: article
license: true
toc: true
key: a20200704
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
julia在计算速度方面要比python快很多，但是在画图方面由于是诞生不久，所以可用的库函数还是比较少，但是可以通过调用python的库函数进行绘图，取长补短。
<!--more-->
# 安装流程
- 安装julia
> 这一步请自行解决，Windows就是点点鼠标，确定自己软件的安装目。Linux就更方便了
- 安装anaconda
> anaconda的安装也自行百度解决，安装它的目的是为了使用jupyter来使用julia和python，这样极其方便(当然了你也可以不赞同，用你喜欢的方式)

安装好anaconda后会连同python一起安装，这是绑定的没所以不用再去考虑单独安装python，以我的安装为例，你可以在目录下看到python
![png](/assets/images/20200704/python.png)
记住这个python的路径，比如我的**D:\\anaconda3\\python.exe**，下面要设置julia的环境变量，使它能正确的找到你的python在哪里，从而可以正确的调用python来作图
```julia
ENV["PYTHON"] = "D:\\anaconda3\\python.exe"
# 这里记得要用双斜杠进行转义操作，不然会提示出错
```
设置好之后就可以进行PyPlot的安装了，至于julia的包如何安装，这一点百度即可，这里不进行赘述。但是需要提醒一些问题，这也是我在安装时候所遇到的。如果安装的时候出现错误，那么首先要保证你已经设置好了上面的python环境，其次的画你要保证你使用anaconda安装了python的作图包[Matplotlib](https://matplotlib.org/)。
anaconda中给python安装包的方式如下，首先在系统中找到你anaconda的命令行，然后进入命令行执行**pip install package-name**，这样你就可以成功的为anaconda中的python安装外部库了。
![png](/assets/images/20200704/anaconda.png)
在确定了python的绘图库以及julia的环境变量中python的位置也设置成功的画，可以进行PyPlot的安装，在Julia中执行**pkg.add("PyPlot")**，在使用Pkg的时候，记得先导入(import Pkg)，之后就可以进行安装了，不过可能还会遇到问题，不过这时候可以执行**Pkg.build("PyCall")**，重新编译一下调用python的这个包，然后再安装PyPlot，应该就可以完整安装了，想在程序中使用的话，调用**using PyPlot**即可。如果是再using过程中出现错误的话，同样可以使用上面提到的重新编译库函数的方式来修复错误，这个信息你也可能会在错误信息提示中看到，可以自己摸索。差不多上面就是我在安装过程中遇到和解决问题的方式。