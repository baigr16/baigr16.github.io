---
title: $\LaTeX$中插入代码
tags: Latex
layout: article
license: true
toc: true
key: a20221215
pageview: true
# cover: /assets/images/GroupTheory/cube_symmetry.jpg
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
最近在适用自己的$\LaTeX$模版整理笔记的时候，想要把自己的笔记和代码整理到一起，这样就不会时间久了笔记找得到反而代码找不到了，整理到一起就更方便了，这里就找到了一个在$\LaTeX$中插入代码的包。
{:.info}
<!--more-->

```latex
\documentclass[11pt]{article}  
\usepackage{pythonhighlight}

\begin{document}
这里插入python代码
\begin{python}
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
config = {
"font.size": 30,
"mathtext.fontset":'stix',
"font.serif": ['SimSun'],
}
rcParams.update(config) # Latex 字体设置
#---------------------------------------------------------
def scatterplot1(cont):
    dataname = "band-theta=" + format(cont,'.2f') + ".dat"
    # dataname = "band-theta=1.05.dat"
    tit = "Twist angele is " + format(cont,'.2f') + "$^\degree$"
    # da1 = "did-short.dat"
    picname = os.path.splitext(dataname)[0] + ".png"
    os.chdir(os.getcwd())# 确定用户执行路径
    x0 = np.loadtxt(dataname)
    plt.figure(figsize=(10,10))
    plt.plot(x0[:,0], x0[:,1:-1], c = 'blue')
    x0min = np.min(x0[:,0])
    x0max = np.max(x0[:,0])
    # y0min = np.min(x0[:,1])
    # y0max = np.max(x0[:,1])
    font2 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 25,
             }
    plt.xlim(x0min,x0max)
    plt.ylim(-2,2)
    # plt.xlabel(r'$\phi_{R-L}/\pi$',font2)
    plt.ylabel("$E$",font2)
    plt.title(tit,font2)
    # plt.yticks(fontproperties='Times New Roman', size = 15)
    plt.xticks(fontproperties='Times New Roman', size = 15)
    # plt.xticks([0,1,2],fontproperties='Times New Roman', size = 25)
    # plt.yticks([-0.25,-0.15,0,0.15,0.25],fontproperties='Times New Roman', size = 25)
    plt.yticks([-1,0,1],fontproperties='Times New Roman', size = 25)
    plt.tick_params(axis='x',width = 2,length = 10)
    plt.tick_params(axis='y',width = 2,length = 10)
    ax = plt.gca()
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5) 
    ax.spines["right"].set_linewidth(1.5)
    ax.spines["top"].set_linewidth(1.5)
    plt.savefig(picname, dpi = 100, bbox_inches = 'tight',transparent = True)
    plt.close()
#---------------------------------------------------------
if __name__=="__main__":
    # scatterplot1(1.05)
    # scatterplot1(0.5)
    scatterplot1(5)
\end{python}
\end{document}
```

虽然这里的环境块为`\python`，在我实际的测试中插入`Julia`的代码也是没有问题的。



