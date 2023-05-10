---
title: WannierTools安装
tags:  vasp
layout: article
license: true
toc: true
key: a20210103b
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
在这里详细的介绍一下如何在Linux中安装WannierTools,这里我要默认服务器上已经安装好了Intel Fortran
{:.info}
<!--more-->
# 安装过程
首先先去[WannierTools](https://www.wanniertools.org/)的官网上面下载好源代码,它的打包格式时`.zip`的,可以使用`unzip filename.zip`来对压缩文件进行解压,如下图所示

![png](/assets/images/wannierTools/s1.png)

解压之后进入到解压的文件目录中,主要关注两个地方,第一个时`INSTALL`这个文件,里面教你怎么安装,第二个就是`src`这个文件夹,里面就是程序包的主要源代码,也包括了编译需要的`Makefile`文件,如下图所示

![png](/assets/images/wannierTools/s2.png)

![png](/assets/images/wannierTools/s3.png)

接下来就是按照`INSTALL`文件中的内容来进行安装,这里我的服务器上已经安装好了并行版本的`Intel Fortran`,所以就选择了并行编译的这个`Makefile`,然后进行编译

![png](/assets/images/wannierTools/s4.png)

经过一段时间的等待之后,就可以在`src`这个文件夹中得到最后的可执行文件`wt.x`,按照文档的指示将`wt.x`复制到`bin`这个文件夹中,到底所有的编译任务就全部完成了,还是很简单的

![png](/assets/images/wannierTools/s5.png)

接下来就是最后一步,将这个可执行文件的路径追加到你的用户`PATH`中,这是`Linux`的内容,想了解的话可以自行百度,这里只要跟着做应该是不会有什么问题

![png](/assets/images/wannierTools/s6.png)

如上图所示,第一个时最终编译得到`wt.x`放置的路径,我们只需要在自己的`.bashrc`这个文件中,执行

```shell
export PATH=$PATH:/home/yxli/opt/wannier_tools-master/bin
```

将这个执行文件的路径追加到`PATH`中即可,然后就可以完整的使用这个工具了.关于这个工具的使用手册和一些技术细节,可以参考`doc`这个文件夹中的内容.以上就是`wannierTools`安装的全部过程了.

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)