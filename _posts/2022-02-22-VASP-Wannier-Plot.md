---
title: 绘制VASP能带和Wannier插值能带进行比较
tags: Topology vasp Python
layout: article
license: true
toc: true
key: a20220222
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
 通常利用Wannier对VASP计算的能带进行拟合的时候，最终总是要看一下到底拟合程度有多好，这里就用Python写个脚本，将二者放在一起进行比较。
{:.info}
<!--more-->
# 数据准备
首先VASP计算的能带结果可以通过[vaspkit](https://vaspkit.com/)来获得，最终可以得到一个`BAND.dat`的数据，使用Wannier90最终也可以得到一个`wannier90_band.dat`的数据，
有了这两个文件之后，就可以进行能带绘制的，下面就是代码。

# 代码
```python
# 通过读取Wannier90给出的能带数据wannier90_band.dat，与VASP结合vaspkit得到的数据BAND.dat
# 来绘制 wannier插值后的能带与VASP计算的能带之间的差距
#-----------------------------------------------------------------
from cProfile import label
import os
import numpy as np
import matplotlib.pyplot as plt
#-----------------------------------------------------------------
def dataread():
    da1 = "wannier90_band.dat"
    x0 = []
    y0 = []
    with open(da1) as file:
            da = file.readlines()
            for f1 in da:
                temp = f1.strip().split() # 移除字符串头尾指定的字符（默认为空格）,并以空格将这些字符串分开
                if len(temp) != 0:
                    x0.append(float(temp[0]))
                    y0.append(float(temp[1]))
    # plt.scatter(x0,y0,s = 1, color = 'red', alpha = 0.7,marker = '.')
    plt.plot(x0,y0,'ro',markersize = 1,label="Wannier90")

    da2 = "BAND.dat"
    fermi =  4.8774
    x1 = []
    y1 = []
    with open(da2) as file:
            da = file.readlines()
            for f1 in da:
                temp = f1.strip().split() # 移除字符串头尾指定的字符（默认为空格）,并以空格将这些字符串分开
                if len(temp) == 2:
                    x1.append(float(temp[0]))
                    y1.append(float(temp[1]) + fermi)
    # plt.scatter(x1,y1 ,s = 1, color = 'blue', alpha = 0.7, marker = '.')
    plt.plot(x1,y1,'bo',markersize = 1,label = "VASP")
    plt.yticks(fontproperties='Times New Roman', size = 15)
    plt.xticks(fontproperties='Times New Roman', size = 15)
    plt.xlim(np.min(x0),np.max(x0))
    plt.legend(loc = 0,ncol = 3)
    plt.ylim(-5,5)
    picname = "band-compare.png"
    plt.savefig(picname, dpi = 600, bbox_inches = 'tight')
#--------------------------------------------------------------
def main():
    os.chdir(os.getcwd())# 确定用户执行路径
    dataread()
#----------------------------------------------------------------
if __name__=='__main__':
    main()
```
![png](/assets/images/python/fig2.png)

这里出了用python之外，还可以使用`gnuplot`来绘制结果

```shell
set style data dots
set encoding iso_8859_1
set terminal png truecolor enhanced font ",50" size 1920, 1680
set size 0.9,1
set output 'band.png'
#set nokey
set key font "Times,24,Bold"
#set key Right at 5,5
set key Right at 3.5,4
#set xrange [0: 4.27065]
#set yrange [-5 : 5]
#set arrow from  0.84341,  -6.29553 to  0.84341,  14.22923 nohead
#set arrow from  2.46594,  -6.29553 to  2.46594,  14.22923 nohead
#set arrow from  3.42724,  -6.29553 to  3.42724,  14.22923 nohead
#set ytics offset -1,0 font 'Times,Bold'
#set xtics offset 0,-0.1 font 'Times,Bold'
set ytics font 'Times,Bold'
set xtics font 'Times,Bold'
#set xtics ("G"  0.00000,"Z"  0.84341,"F"  2.46594,"G"  3.42724,"L"  4.27065)
#plot "wannier90_band.dat" w p  pt 7  ps 1.1 lc 'red',"BAND.dat" w p pt 7 ps 1.1 lc 'blue'
fermi =  4.8774
plot "wannier90_band.dat" w l  lw 5.0 lt 7 lc 'red' title "wannier","BAND.dat" u 1:($2+fermi) w l lw 5.0 lt 7  lc 'blue' title "VASP"
```

![png](/assets/images/python/fig3.png)