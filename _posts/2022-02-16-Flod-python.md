---
title: Python实现文件夹中文件批量处理
tags: Python Study
layout: article
license: true
toc: true
key: a20220216
pageview: true
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
Python实现文件夹中问价批量复制并移动到其他位置。
{:.info}
<!--more-->

# 问题描述
一个文件中有许多个子文件夹，而子文件夹下面又有名称相同的一个文件夹，形成了一层多余的结果，所以想把最下面的文件提取出来，只保留最外层的这个文件夹下面就可以，实现代码如下，这里比较懒，没有进行函数化
```python
from distutils import filelist
import os
import shutil
#--------------------------------------------------------
os.chdir(os.getcwd())# 确定用户执行路径
tarfolename = "newfold" # 新建一个文件夹来存储这些文件，消除掉原本有两层的结构
folder = os.path.exists(tarfolename)
if not folder:
    os.mkdir(tarfolename) # 在当前路径下新建一个文件夹，用来存放新的文件
else:
    print("文件夹已经存在")

filelist = os.listdir(os.getcwd()) # 获取当前文件夹中的所有文件名称(包括文件夹)

flist2 = []
for filename in filelist:
    f1 = filename  + "/" + filename
    f2 = os.getcwd() + "/" +  f1 # 文件夹下面还有名称相同的文件夹，这里进行拼接，组合成正确的文件路径，方便后面获取其中的文件
    flist2.append(f2)


f11 = os.getcwd() + "/" + tarfolename # 文件复制的目标路径

print(f11)
for f1 in flist2: #开始遍历每个文件夹，这里已经考虑了文件夹下面还有文件夹
    f2 = os.listdir(f1) # 获取文件夹中的所有文件
    for f3 in f2: # 遍历文件夹中的文件
        src = os.path.join(f1,f3) # 文件路径

        folder = os.path.exists(f11 + "/" + os.path.basename(f1)) # 先判断文件夹是否存在

        if not folder: # 如果不存在就新建这个文件夹，然后才能copy文件进来
            os.mkdir(f11 + "/" + os.path.basename(f1)) # 在当前路径下新建一个文件夹，用来存放新的文件
        else:
            print("文件夹已经存在")

        tar = os.path.join(f11 + "/" + os.path.basename(f1),f3) # 需要移动到哪个位置

        shutil.move(src,tar)
```
![png](/assets/images/python/fig1.jpg)

# 结果

最终的效果就是把本来套了两层文件夹中的数据，去掉了中间的一层，只保留了最外层的文件夹名称。