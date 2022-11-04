---
title: Jupyter中安装Julia
tags: Study Fortran
layout: article
license: true
toc: true
pageview: true
key: a202008061
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
Jupyter在数据科学方面用起来还是很方便，而且在里面可以同时运行多种变成语言，正好最近换了电脑需要重新安装julia，也就正好记录一下安装过程。
{:.success}
<!--more-->
Jupyter是在Ananconda的组件，所以只要先下载安装[Ananconda](https://www.anaconda.com/)，那么自然就可以看到jupyter，且在安装过程中会自动安装好python，就不用再单独安装python了，而且此时的jupyter中是自带python的，接下来就是安装julia了。

[julia](https://cn.julialang.org/downloads/)在中文社区下载好之后，直接运行安装即可，这里记一下自己的安装目录，比如我的安装目录是**D:\Julia-1.4.1**，安装完成后就可以正常启动julia了，不过此时的julia只不过是一个基本的内容，还需要安装一系列的包。

在安装包之前要说明一下，我们在这里需要使用一下镜像网站来安装包，否则速度会非常慢，而且经常会失败，这里用到的是北京外国语的一个镜像，首先要找到你julia的**startup.jl**文件，这个文件的路径为**D:\Julia-1.4.1\etc\julia\startup.jl**也就是在安装目录下寻找，找到之后，刚开始这个文件中只有一些注释过的语句，我们需要在这个文件中加入
> **ENV["JULIA_PKG_SERVER"] = "https://mirrors.bfsu.edu.cn/julia/static"**

加入之后，重启julia，然后查看版本信息即如图所示
![png](/assets/images/Julia/julia-jupyter.png)
看到上图所示的信息之后就说明镜像设置成功了，接下就是安装将julia嵌入Jupyter的包IJulia了
> import Pkg
Pkg.add("IJulia")
using IJulia
 
以此执行上面的三个命令之后，即可成功在Jupyter中安装Julia，接下来打开Jupyter进行新建文件是，就可以看到Julia的命令了
![png](/assets/images/Julia/julia2.png)