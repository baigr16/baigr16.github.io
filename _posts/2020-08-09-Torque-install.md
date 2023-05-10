---
title: 队列系统Torque安装
tags: Linux Shell
layout: article
license: true
toc: true
key: a202008091
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
在进行第一性原理计算的时候，是需要用到队列系统去提交任务的，在这里就记录一下自己安装队列系统的过程
<!--more-->
# 安装
一般情况下Linux的软件都安装到/opt这个文件夹中，所以先在这里新建文件夹
> cd /opt    mkdir torque     cd torque

首先你要确定自己有安装包，我这里的安装包是_Torque-6.1.2_，先将软件进行解压`tar zxvf torque-6.1.2.tar.gz`，接下来执行`./configure --with-default-server=master`，这里默认都是利用root进行的。这个设置完成后会提示提可以进行make，那么就执行`make`。如果没有安装包可执行
> wget http://wpfilebase.s3.amazonaws.com/torque/torque-6.1.2.tar.gz
tar -zxvf torque-6.1.2.tar.gz
cd torque-6.1.2/

安装一些必要的依赖
> yum install libxml2-devel openssl-devel gcc gcc-c++ boost-devel libtool-y

依赖安装完成后，进入到torque-6.1.2中执行
> ./configure --prefix=/usr/local/torque --with-scp--with-default-server=master
make
make install
make packages

这三个命令都需要执行一定的时间。接下来执行
> libtool --finish /usr/local/torque/lib
## 服务配置
> cp contrib/init.d/{pbs_{server,sched,mom},trqauthd} /etc/init.d/
for i in pbs_server pbs_sched pbs_mom trqauthd; do chkconfig --add $i; chkconfig $i on; done
## Torque环境变量设置
> TORQUE=/usr/local/torque
echo "TORQUE=$TORQUE" >> /etc/profile
echo "export PATH=\$PATH:$TORQUE/bin:$TORQUE/sbin" >> /etc/profile
source /etc/profile
# Ref
[Blog](https://blog.csdn.net/zhaosongbin/article/details/87914746)

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)