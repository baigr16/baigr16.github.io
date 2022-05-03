---
title: Kwant学习
tags: transport Python
layout: article
license: true
toc: true
key: a20220503
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
最近打算转学一下输运,首先项借助`Kwant`学习一下其中的例子先帮助自己理解一些输运中的概念和基本的计算,后面再自己利用格林函数的方法边学习边实践,这里还是想整理记录一下自己在学习`Kwant`的实例中的一些代码注释和笔记.
{:.info}
<!--more-->
# 正方晶格散射
```python
from matplotlib import pyplot
import kwant

# 首先定义一个系统

syst = kwant.Builder()


a = 1 
lat = kwant.lattice.square(a) # 确定系统的格点位置,这里选择的是正方格点,那么此时lat就是一个正方点阵

t = 1.0
W = 10
L = 30

# 通过实空间的哈密顿量来确定散射区域

for i in range(L):
    for j in range(W):
        # On-site Hamiltonian
        syst[lat(i, j)] = 4 * t # 通过格点来确定hopping以及onsite占位能

        # Hopping in y-direction
        if j > 0:
            syst[lat(i, j), lat(i, j - 1)] = -t

        # Hopping in x-direction
        if i > 0:
            syst[lat(i, j), lat(i - 1, j)] = -t # 通过这个赋值可以看到在边界上是不存在hopping的,所以这些位置就是hard-wall

# Then, define and attach the leads:

# 这里定义一个电极,通过矢量
sym_left_lead = kwant.TranslationalSymmetry((a, 0)) # 通过矢量(a,0)来定义,此时这个方向一定要是远离散射区域的方向,那么也就表示着从散射区域的右端都是电极
left_lead = kwant.Builder(sym_left_lead) # 定义一个电极

for j in range(W):
    left_lead[lat(0, j)] = 4 * t # 这个相当于是电极上和散射区域具有相同的占位能
    if j > 0: 
        left_lead[lat(0, j), lat(0, j - 1)] = -t # 电极与散射区域的边界上的耦合
    left_lead[lat(1, j), lat(0, j)] = -t

syst.attach_lead(left_lead) # 将电极和散射区域连接起来

# 安装完了左边的电极再来定义并安装右边的电极
sym_right_lead = kwant.TranslationalSymmetry((-a,0))
right_lead = kwant.Builder(sym_right_lead)

for j in range(W):
    right_lead[lat(0, j)] = 4 * t
    if j > 0:
        right_lead[lat(0, j), lat(0, j - 1)] = -t
    right_lead[lat(1, j), lat(0, j)] = -t

syst.attach_lead(right_lead) # 把右边的电极也安装到系统上


kwant.plot(syst) # 看一下构建的系统长什么模样


syst = syst.finalized() # 初始化系统,后面就可以对系统的各种性质进行计算了

# Now that we have the system, we can compute conductance
energies = []
data = []
for ie in range(100):
    energy = ie * 0.01

    # 在确定能量下面计算散射矩阵
    smatrix = kwant.smatrix(syst, energy)

    # compute the transmission probability from lead 0 to
    # lead 1
    energies.append(energy)
    # 在前面我们分别在左右两端都定义了电极,它们都是由相对应的编号的,此时就可以通过编号来计算从0电极到1电极上的投射率
    # 熟悉输运的话就知道我们可以通过散射矩阵计算得到投射率
    data.append(smatrix.transmission(1, 0)) 

# 画图展示电导
pyplot.figure()
pyplot.plot(energies, data)
pyplot.xlabel("energy [t]")
pyplot.ylabel("conductance [e^2/h]")
pyplot.show()
```
上面的代码详细的将构造的细节都展现了出来,但是同样发现有些代码是重复的,这里借助`Kwant`的一些内置属性以及`python`的一些技巧将代码简化
```python
import kwant

# For plotting
from matplotlib import pyplot


def make_system(a=1, t=1.0, W=10, L=30):
    # Start with an empty tight-binding system and a single square lattice.
    # `a` is the lattice constant (by default set to 1 for simplicity.
    lat = kwant.lattice.square(a)

    syst = kwant.Builder()
    syst[(lat(x, y) for x in range(L) for y in range(W))] = 4 * t
    syst[lat.neighbors()] = -t
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))
    lead[(lat(0, j) for j in range(W))] = 4 * t
    lead[lat.neighbors()] = -t # 通过Kwant内建的属性得到最近邻的hopping,这里neighbors()没有给值就是默认最近邻,然也可以给值来确定是第几近邻
    # 定义电极
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))
    lead[(lat(0, j) for j in range(W))] = 4 * t
    lead[lat.neighbors()] = -t
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed()) #因为左右电极只是电极平移不变的方向不一样,所以这里直接使用了reversed这个属性

    return syst
#---------------------------------
def plot_conductance(syst, energies):
    # 这就是个画图程序
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()
#-----------------------------------
def main():
    syst = make_system()

    # Check that the system looks as intended.
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # We should see conductance steps.
    plot_conductance(syst, energies=[0.01 * i for i in range(100)])
#-------------------------------------
if __name__ == '__main__':
    main()
```

# 考虑自旋
```python
import kwant

# For plotting
from matplotlib import pyplot

# For matrix support,这个包在小矩阵运算方面速度较快
import tinyarray

# define Pauli-matrices for convenience
sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])


def make_system(t=1.0, alpha=0.5, e_z=0.08, W=10, L=30):
    # 定义方格点,这里没有给定晶格长度,但是系统默认为1
    lat = kwant.lattice.square()

    syst = kwant.Builder() # 构建散射区域对象

    # 这里在定义散射区域的时候,每个位置上的hopping大小同样可以利用矩阵表示
    syst[(lat(x, y) for x in range(L) for y in range(W))] = \
        4 * t * sigma_0 + e_z * sigma_z
    # hoppings in x-direction
    syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_y / 2 # x正方向hopping
    # hoppings in y-directions
    syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_x / 2 # y正方向hopping,负方向会默认设置

    #定义左端电极
    lead = kwant.Builder(kwant.TranslationalSymmetry((-1, 0)))

    lead[(lat(0, j) for j in range(W))] = 4 * t * sigma_0 + e_z * sigma_z 
    # hoppings in x-direction
    lead[kwant.builder.HoppingKind((1, 0), lat, lat)] = \
        -t * sigma_0 + 1j * alpha * sigma_y / 2
    # hoppings in y-directions
    lead[kwant.builder.HoppingKind((0, 1), lat, lat)] = \
        -t * sigma_0 - 1j * alpha * sigma_x / 2

    # 将电极和散射区域连接起来
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst


def plot_conductance(syst, energies):
    # Compute conductance
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()


def main():
    syst = make_system()

    # 看一下设置的系统样式
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # 画图
    plot_conductance(syst, energies=[0.01 * i - 0.3 for i in range(100)])


if __name__ == '__main__':
    main()
```

# 空间变化势能
```python
import kwant

# For plotting
from matplotlib import pyplot


def make_system(a=1, t=1.0, W=10, L=30, L_well=10):
    #构建正方点阵,晶格参数为a
    lat = kwant.lattice.square(a)

    syst = kwant.Builder() # 构建系统

    #定义散射区域,此时散射区域中的势能是随空间位置变化的

    def potential(site, pot):
        (x, y) = site.pos
        if (L - L_well) / 2 < x < (L + L_well) / 2:
            return pot
        else:
            return 0

    def onsite(site, pot):
        return 4 * t + potential(site, pot)

    syst[(lat(x, y) for x in range(L) for y in range(W))] = onsite
    syst[lat.neighbors()] = -t

    #定义电极并将其与散射区域耦合
    lead = kwant.Builder(kwant.TranslationalSymmetry((-a, 0)))
    lead[(lat(0, j) for j in range(W))] = 4 * t
    lead[lat.neighbors()] = -t
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst


def plot_conductance(syst, energy, welldepths):

    # Compute conductance
    data = []
    for welldepth in welldepths:
        smatrix = kwant.smatrix(syst, energy, params=dict(pot=-welldepth))
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(welldepths, data)
    pyplot.xlabel("well depth [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()


def main():
    syst = make_system()

    # 查看定义的系统
    kwant.plot(syst)
    syst = syst.finalized()

    # 计算电导
    plot_conductance(syst, energy=0.2,
                     welldepths=[0.01 * i for i in range(100)])

if __name__ == '__main__':
    main()
```

# 加入磁通
```python
from cmath import exp
from math import pi

import kwant

# For plotting
from matplotlib import pyplot


def make_system(a=1, t=1.0, W=10, r1=10, r2=20):
    lat = kwant.lattice.square(a) # 定义四方点阵

    syst = kwant.Builder() # 构建系统

    #定义散射区域
    # 这里想定义一个圆环区域
    def ring(pos):
        (x, y) = pos
        rsq = x ** 2 + y ** 2
        return (r1 ** 2 < rsq < r2 ** 2)

    # 通过一个函数来限制格点的位置
    syst[lat.shape(ring, (0, r1 + 1))] = 4 * t
    syst[lat.neighbors()] = -t

    # 这里想要加磁场,所以每个hopping上面就要加入磁场带来的位相
    def hopping_phase(site1, site2, phi):
        return -t * exp(1j * phi)

    def crosses_branchcut(hop):
        ix0, iy0 = hop[0].tag

        # builder.HoppingKind with the argument (1, 0) below
        # returns hoppings ordered as ((i+1, j), (i, j))
        return iy0 < 0 and ix0 == 1  # ix1 == 0 then implied

    # 此时要考虑磁场,就需要对hopping大小重新进行修改,选择确定的磁场规范之后即只修改一个方向的hopping就可以
    def hops_across_cut(syst):
        for hop in kwant.builder.HoppingKind((1, 0), lat, lat)(syst):
            if crosses_branchcut(hop):
                yield hop
    syst[hops_across_cut] = hopping_phase

    #### Define the leads. ####
    # left lead
    sym_lead = kwant.TranslationalSymmetry((-a, 0))
    lead = kwant.Builder(sym_lead)

    def lead_shape(pos):
        (x, y) = pos
        return (-W / 2 < y < W / 2)

    lead[lat.shape(lead_shape, (0, 0))] = 4 * t
    lead[lat.neighbors()] = -t

    #### Attach the leads and return the system. ####
    syst.attach_lead(lead)
    syst.attach_lead(lead.reversed())

    return syst


def plot_conductance(syst, energy, fluxes):
    # compute conductance

    normalized_fluxes = [flux / (2 * pi) for flux in fluxes]
    data = []
    for flux in fluxes:
        smatrix = kwant.smatrix(syst, energy, params = dict(phi = flux))
        data.append(smatrix.transmission(1, 0))

    pyplot.figure()
    pyplot.plot(normalized_fluxes, data)
    pyplot.xlabel("flux [flux quantum]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()


def main():
    syst = make_system()

    # 查看定义的系统的形状
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    # W这里通过改变磁通大小来计算电导变化
    plot_conductance(syst, energy=0.15, fluxes=[0.01 * i * 3 * 2 * pi
                                                for i in range(100)])

if __name__ == '__main__':
    main()
```

# 电极能带模式计算
```python
import kwant
from matplotlib import pyplot

def make_lead(a=1, t=1.0, W=10): # 定义一个电极
    
    lat = kwant.lattice.square(a)

    sym_lead = kwant.TranslationalSymmetry((-a, 0))
    lead = kwant.Builder(sym_lead)

    # build up one unit cell of the lead, and add the hoppings
    # to the next unit cell
    for j in range(W):
        lead[lat(0, j)] = 4 * t

        if j > 0:
            lead[lat(0, j), lat(0, j - 1)] = -t

        lead[lat(1, j), lat(0, j)] = -t # 元胞间hopping

    return lead


def main():
    lead = make_lead().finalized()
    kwant.plotter.bands(lead, show=False)
    pyplot.xlabel("momentum [(lattice constant)^-1]")
    pyplot.ylabel("energy [t]")
    pyplot.show()

if __name__ == '__main__':
    main()
```
# 闭合系统电流计算
```python
from cmath import exp
import numpy as np
import kwant
from matplotlib import pyplot
import scipy.sparse.linalg as sla # 稀疏矩阵本征值计算



def make_system(a=1, t=1.0, r=10):
    lat = kwant.lattice.square(a, norbs=1)

    syst = kwant.Builder()

    # 定义一个圆形的量子点区域
    def circle(pos):
        (x, y) = pos
        rsq = x ** 2 + y ** 2
        return rsq < r ** 2

    def hopx(site1, site2, B):
        # 这里在定义hopping的时候有一个参数B,后面在构建系统的时候就会通过字典的形式给值
        y = site1.pos[1]
        return -t * exp(-1j * B * y)

    syst[lat.shape(circle, (0, 0))] = 4 * t
    # x方向hopping
    syst[kwant.builder.HoppingKind((1, 0), lat, lat)] = hopx
    # y方向hopping
    syst[kwant.builder.HoppingKind((0, 1), lat, lat)] = -t
    # 没有设计电极,所以这是一个闭合的系统
    return syst


def plot_spectrum(syst, Bfields):

    energies = []
    for B in Bfields:
        # 这里就会通过参数B给定确定的值
        ham_mat = syst.hamiltonian_submatrix(params=dict(B=B), sparse=True)

        # 计算最低的15条本征值,这是对角化函数的设置
        ev = sla.eigsh(ham_mat.tocsc(), k = 15, sigma = 0,return_eigenvectors = False)
        energies.append(ev)

    pyplot.figure()
    pyplot.plot(Bfields, energies)
    pyplot.xlabel("magnetic field [arbitrary units]")
    pyplot.ylabel("energy [t]")
    pyplot.show()

def sorted_eigs(ev):
    evals, evecs = ev
    evals, evecs = map(np.array, zip(*sorted(zip(evals, evecs.transpose()))))
    return evals, evecs.transpose()

def plot_wave_function(syst, B=0.001):
    # 计算系统的波函数
    ham_mat = syst.hamiltonian_submatrix(sparse=True, params=dict(B=B))
    evals, evecs = sorted_eigs(sla.eigsh(ham_mat.tocsc(), k=20, sigma=0)) # k=20表示计算能量最低的20个模式

    # 绘制波函数几率分布
    kwant.plotter.map(syst, np.abs(evecs[:, 9])**2,
                      colorbar=True, oversampling=1)


def plot_current(syst, B=0.001):
    # Calculate the wave functions in the system.
    ham_mat = syst.hamiltonian_submatrix(sparse=True, params=dict(B=B))
    evals, evecs = sorted_eigs(sla.eigsh(ham_mat.tocsc(), k=20, sigma=0))

    # 计算电流
    J = kwant.operator.Current(syst) 
    current = J(evecs[:, 9], params=dict(B=B))
    kwant.plotter.current(syst, current, colorbar=True)


def main():
    syst = make_system()

    # Check that the system looks as intended.
    kwant.plot(syst)

    # Finalize the system.
    syst = syst.finalized()

    try:
        plot_spectrum(syst, [iB * 0.002 for iB in range(100)])

        syst = make_system(r=30).finalized()
        plot_wave_function(syst)
        plot_current(syst)
    except ValueError as e:
        if e.message == "Input matrix is not real-valued.":
            print("The calculation of eigenvalues failed because of a bug in SciPy 0.9.")
            print("Please upgrade to a newer version of SciPy.")
        else:
            raise

if __name__ == '__main__':
    main()
```

# 石墨烯结构计算
```python
from math import pi, sqrt, tanh
import kwant
import scipy.sparse.linalg as sla
from matplotlib import pyplot



sin_30, cos_30 = (1 / 2, sqrt(3) / 2)
graphene = kwant.lattice.general([(1, 0), (sin_30, cos_30)],
                                 [(0, 0), (0, 1 / sqrt(3))]) # 定义石墨烯的格点矢量
a, b = graphene.sublattices;


def make_system(r=10, w=2.0, pot=0.1):

    # 定义散射区域,这是个圆形区域
    def circle(pos):
        x, y = pos
        return x ** 2 + y ** 2 < r ** 2

    syst = kwant.Builder() # 构建系统

    # 空间变化的是能
    def potential(site):
        (x, y) = site.pos
        d = y * cos_30 + x * sin_30
        return pot * tanh(d / w)

    syst[graphene.shape(circle, (0, 0))] = potential  # 对给定形状的散射区域设置势能

    # 设置hopping
    hoppings = (((0, 0), a, b), ((0, 1), a, b), ((-1, 1), a, b))
    syst[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    # 将前面通过hopping设置的区域进行修改,这里是删除掉了一些位置,并对某几个位置的hopping进行了修改
    del syst[a(0, 0)]
    syst[a(-2, 1), b(2, 2)] = -1

    #定义电极
    # left lead
    sym0 = kwant.TranslationalSymmetry(graphene.vec((-1, 0)))

    def lead0_shape(pos):
        x, y = pos
        return (-0.4 * r < y < 0.4 * r)

    lead0 = kwant.Builder(sym0)
    lead0[graphene.shape(lead0_shape, (0, 0))] = -pot
    lead0[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    # 定义右电极
    sym1 = kwant.TranslationalSymmetry(graphene.vec((0, 1)))

    def lead1_shape(pos):
        v = pos[1] * sin_30 - pos[0] * cos_30
        return (-0.4 * r < v < 0.4 * r)

    lead1 = kwant.Builder(sym1)
    lead1[graphene.shape(lead1_shape, (0, 0))] = pot
    lead1[[kwant.builder.HoppingKind(*hopping) for hopping in hoppings]] = -1

    return syst, [lead0, lead1]


def compute_evs(syst):
    # Compute some eigenvalues of the closed system
    sparse_mat = syst.hamiltonian_submatrix(sparse=True)

    evs = sla.eigs(sparse_mat, 2)[0]
    print(evs.real)


def plot_conductance(syst, energies):
    # 计算从0端到1端的电导
    data = []
    for energy in energies:
        smatrix = kwant.smatrix(syst, energy)
        data.append(smatrix.transmission(0, 1))

    pyplot.figure()
    pyplot.plot(energies, data)
    pyplot.xlabel("energy [t]")
    pyplot.ylabel("conductance [e^2/h]")
    pyplot.show()


def plot_bandstructure(flead, momenta):
    bands = kwant.physics.Bands(flead)
    energies = [bands(k) for k in momenta]

    pyplot.figure()
    pyplot.plot(momenta, energies)
    pyplot.xlabel("momentum [(lattice constant)^-1]")
    pyplot.ylabel("energy [t]")
    pyplot.show()


def main():
    pot = 0.1
    syst, leads = make_system(pot=pot)

    # To highlight the two sublattices of graphene, we plot one with
    # a filled, and the other one with an open circle:
    def family_colors(site):
        return 0 if site.family == a else 1

    # Plot the closed system without leads.
    kwant.plot(syst, site_color=family_colors, site_lw=0.1, colorbar=False)

    # Compute some eigenvalues.
    compute_evs(syst.finalized())

    # Attach the leads to the system.
    for lead in leads:
        syst.attach_lead(lead)

    # Then, plot the system with leads.
    kwant.plot(syst, site_color=family_colors, site_lw=0.1,
               lead_site_lw=0, colorbar=False)

    # Finalize the system.
    syst = syst.finalized()

    # Compute the band structure of lead 0.
    momenta = [-pi + 0.02 * pi * i for i in range(101)]
    plot_bandstructure(syst.leads[0], momenta)

    # Plot conductance.
    energies = [-2 * pot + 4. / 50. * pot * i for i in range(51)]
    plot_conductance(syst, energies)


if __name__ == '__main__':
    main()
```