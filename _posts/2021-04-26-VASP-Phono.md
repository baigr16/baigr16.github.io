---
title: VASP+Phonopy计算声子谱
tags: vasp
layout: article
license: true
toc: true
key: a20210426
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
这篇博客是我学习声子谱计算的一些笔记,因为通常在判断一个体系是否具有稳定性的时候,需要计算其声子谱,最近也在慢慢摸索第一性计算的相关工具和知识,就一起整理出来.
<!--more-->

#  Phonopy 安装
VASP的安装这里就不多说了,可以参考[VASP编译安装](https://yxli8023.github.io/2020/08/09/VASP-install.html)这篇博客中的内容,这里主要先整理如何安装[Phonopy](https://phonopy.github.io/phonopy/), 完全参考的是[官网](https://phonopy.github.io/phonopy/)上的教程,从官网上的安装教程来看,最好是先安装好[Ananconda](https://www.anaconda.com/),关于Ananconda的安装可以参考[做数值计算好用的软件及杂项整理](https://yxli8023.github.io/2020/09/16/introduction.html)这篇博客中的内容.在安装好了Ananconda之后,开始安装Phonopy
```shell
conda install -c conda-forge phonopy
```
我是在Linux服务器上安装的,所以在安装Ananconda与Phonopy时都是以root用户进行的,中间只会简单的进行库函数的更新与Phonopy的安装,过程耗时非常少,耐心等待几分钟即可.

#  自洽
首先在计算声子谱的时候,要先保证完成了结构优化,因为我对结构优化还不是很熟悉,所以我在这里整理的只是如何正确的完整声子谱计算的整个流程,我是利用一个已知的结构,在自洽计算的结果上进行声子谱计算的.
{:.warning}
下面是我做自洽计算的VASP输入文件
- INCAR

```shell
Global Parameters
ISTART =  1            (Read existing wavefunction; if there)
ISPIN =  2           (Spin polarised DFT)
ICHARG =  2         (Non-self-consistent: GGA/LDA band structures)
LREAL  = .FALSE.       (Projection operators: automatic)
ENCUT  =  400        (Cut-off energy for plane wave basis set, in eV)
PREC   =  Normal       (Precision level)
LWAVE  = .TRUE.        (Write WAVECAR or not)
LCHARG = .TRUE.        (Write CHGCAR or not)
ADDGRID= .TRUE.        (Increase grid; helps GGA convergence)

# GGA = PE

#  LVTOT  = .TRUE.      (Write total electrostatic potential into LOCPOT or not)

#  LVHAR  = .TRUE.      (Write ionic + Hartree electrostatic potential into LOCPOT or not)

#  NELECT =             (No. of electrons: charged cells; be careful)

#  LPLANE = .TRUE.      (Real space distribution; supercells)

#  NPAR   = 4           (Max is no. nodes; don't set for hybrids)

#  NWRITE = 2           (Medium-level output)

#  KPAR   = 2           (Divides k-grid into separate groups)

#  NGX    = 500         (FFT grid mesh density for nice charge/potential plots)

#  NGY    = 500         (FFT grid mesh density for nice charge/potential plots)

#  NGZ    = 500         (FFT grid mesh density for nice charge/potential plots)
 
Static Calculation
ISMEAR =  0            (gaussian smearing method)
SIGMA  =  0.05         (please check the width of the smearing)
LORBIT =  11           (PAW radii for projected DOS)
NEDOS  =  2001         (DOSCAR points)
NELM   =  60           (Max electronic SCF steps)
EDIFF  =  1E-08        (SCF energy convergence; in eV)
NBANDS = 64

# NBANDS = 128
LSORBIT    = .TRUE.    (Activate SOC)
```

- POSCAR
```shell
Bi2Se3
1.0
-2.069  -3.583614  0.000000
 2.069  -3.583614  0.000000
 0.000   2.389075  9.546667
Bi   Se
 2   3
Direct
 0.3990    0.3990    0.6970
 0.6010    0.6010    0.3030
 0     0     0.5
 0.2060    0.2060    0.1180
 0.7940    0.7940    0.8820             
```

- KPOINTS
```shell
Monkhorst Pack
0
G
 12  12  3
 0   0   0
```
完整自洽计算之后,下一步就开始构建超胞计算声子谱
```shell
phonopy -d --dim="4 4 1"
```
这里进行扩胞的时候,要让3个方向都大于$10A$的量级,运行这个命令之后会生成很多个POSCAR-001之类的文件,不过我们使用的暂时只有**SPOSCAR**这个文件
{:.info}

```shell
mv POSCAR POSCAR-unit
mv SPOSCAR POSCAR
```
先保存原来的结构文件,之后将超胞的结构文件复制为**POSCAR**,此时查看可以发现元胞中原子数会变得非常的多,所以在计算的时候需要修改**KPOINTS**文件
- KPOINTS
```shell
Monkhorst Pack
0
G
 1   1   1
 0   0   0
```
这里选取一个$\Gamma$点就可以,太大的话可能无法计算,但这里只取一个$\Gamma$点还有其他的原因,我正在学习中.
{:.warning}

其中的INCAR不需要修改,接下来就可以提交计算了
```shell
mpirun -np 24 vasp_std &
```
此时的计算相比于前一步自洽需要更长的时间.

# 公众号
相关内容均会在公众号进行同步，若对该Blog感兴趣，欢迎关注微信公众号。
{:.info}

![png](/assets/images/qrcode.jpg)