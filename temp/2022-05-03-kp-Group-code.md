---
title: 对称性约束$\mathbf{k}\cdot \mathbf{p}$哈密顿量(算法解析)
tags: Group-Theory Mathematica Code 
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
我在前面先预告了一版[对称性约束$\mathbf{k}\cdot \mathbf{p}$哈密顿量(预告版)](https://yxli8023.github.io/2022/04/30/kp-Group.html),这里就详细的整理一下这个代码的思路和每一行代码的想法.
{:.info}
<!--more-->
# 理论背景
其实对称性约束$\mathbf{k}\cdot\mathbf{p}$哈密顿量的核心公式就是

$$
C_mH(\mathbf{k})C_m^{-1}=H(R_m\mathbf{k})\label{q1}
$$

而接下来我们需要知道的是两件事情,第一就是

$$
R_m\mathbf{k}\rightarrow ?
$$

第二就是$C_m$操作的矩阵表示,这两个就是计算所需要的输入参数.

通过公式(\ref{q1})我们可以明白,当确定了操作的表示之后,那么表示的维度就与哈密顿量的维度是相同的了,因为公式(\ref{q1})其实看起来不单单是个抽象的表达式,你也可以将它看做是矩阵运算,这当然是没问题的.

当知道了上面两个输入参数之后,接下来就是如何表示哈密顿量了.其实也很简单,既然我是要求具体的形式,那么我就先将它展开成一个一般的形式,然后在通过对称性进行约束,也就是公式(\ref{q1})的作用了.那么接下来就将哈密顿量展开

$$
H(k_x,k_y,k_z)=1\times\Gamma_1 + k_x\Gamma_2 + k_y\Gamma_3 + k_z\Gamma_4 + k_xk_y\Gamma_5 + k_xk_z\Gamma_6 + k_yk_z\Gamma_7 + \cdots
$$

其实也就将哈密顿量类似于Taylor级数一样展开,只不过这里有两点不同,第一点是多项式前面的展开系数都是1,因为我们最终总是要通过具体的参数来表示前面的系数,所以这里有没有个常数的系数其实并没有什么影响,我们总可以$3\times c_1\rightarrow c_1$.第二点不同就是这里不但有多项式的展开,同样还有矩阵展开$\Gamma_i$,这里的每一个矩阵的维度都与我们选择的表示的维度是相同的.

在程序中我们想要实现这个哈密顿量展开就需要一点小技巧,首先我们知道哈密度量一定是厄米的,所以我们可以将厄米矩阵分解成一个对称矩阵和反对称矩阵的和

$$
H=A+i\times B\\
H^\dagger=A-i\times B\\
A^T=A\quad B^T=-B
$$

这里的$i$就是虚数单位,$A$是个对称矩阵,$B$就是一个反对称矩阵.这么做的优点就是当我们再考虑反幺正操作(操作中带有复共轭操作)的时候,我们可以直接利用$H^\dagger$而不是对$H$去求共轭,这个在我们想法上看起来好像是等价的,但是在程序实现起来就有点麻烦了,因为如果想要直接对$H$求共轭操作,那么涉及到的是符号计算,这在程序执行过程中就会很麻烦.所以分别定义$H^\dagger$和$H$在程序计算的时候就不用单独执行共轭操作了,这是写这个程序的第一个小技巧.

那么通过对称矩阵和反对称矩阵构建出来哈密顿量之后,如何确定那些参数应该消失呢?还是利用公式

$$
C_mH(\mathbf{k})C_m^{-1}=H(R_m\mathbf{k})
$$

我们可以通过符号计算分别得到等式左右两端的表达式,这个时候对两个结果求差,那么公式(\ref{q1})告诉我们在对称性的约束下它们的差一定要是零,所以我们就需要令程序计算中两者的差为零,这样就可以得到哈密顿量展开中,$\Gamma_i$矩阵中参数之间的关系.

比如某一个$\Gamma_i$矩阵中分别有参数

$$c_1,c_2,c_3,c_4$$

在刚开始的时候它们是相互独立的,但是再通过对称性的限制之后可能就有

$$
c_1=c_2=c_3
$$

这样我们就将原本4个独立的参数约化成了两个独立的参数.这也就是利用对称性约束来得到参数约化的目的.而对于每一个对称操作,我们都可以进行这样的约化过程.但是这里在写程序的时候有个小track,虽然不同的对称操作对哈密顿量约束的时候是相互独立的,但是我们的目的就是要利用这些对称操作来将哈密顿量展开的矩阵形式中的参数约束到最少,所以在依次利用对称操作来约束哈密顿量的时候,比如我们用$C_1$约束完了哈密顿量得到了一些参数之间的约束关系

$$
c_1=c_2=c_3
$$

我们要现将这些关系代入到我们哈密顿量的参数中,这样才算是达到了参数约化的目的,此时相当于我们得到了一个新的哈密顿量,接下来就是将系统的每一个对称操作轮流执行一遍,这样就可以得到在这些对称性约束下面的最简形式的哈密顿量.以上也就是实现对称性约束$\mathbf{k}\cdot\mathbf{p}$哈密顿量的主要思想.

# 代码解析

```julia
ClearAll["`*"]
sd1k = {kx -> kx/2 - \[Sqrt]3 *ky/2, ky -> -Sqrt[3] kx/2 - ky/2, 
  kz -> kz}; C3zk = {kx -> -kx/2 + \[Sqrt]3 *ky/2, 
  ky -> -Sqrt[3] kx/2 - ky/2, kz -> kz}; shk = {kx -> kx, ky -> ky, 
  kz -> -kz};
PTk = {kx -> kx, ky -> ky, kz -> kz};
sd1 = -I*KroneckerProduct[PauliMatrix[0], 
   PauliMatrix[
    1]]; C3z = (-KroneckerProduct[PauliMatrix[0], PauliMatrix[0]] - 
    I*Sqrt[3]*KroneckerProduct[PauliMatrix[3], PauliMatrix[3]])/
  2; sh = KroneckerProduct[PauliMatrix[0], PauliMatrix[3]];
PT = KroneckerProduct[PauliMatrix[1], PauliMatrix[0]];
UOperationRes = {{C3z, C3zk}, {sh, shk}, {sd1, sd1k}};(*整合所有的幺正操作*)
AUOperationRes = {{PT, PTk}};(*整合所有的反幺正操作*)
OperationRes = {UOperationRes, 
   AUOperationRes};(*将所有对称操作整合起来,总体上分为两类,幺正操作和非幺正操作*)
(*----------------------------------------------------------------------------------------*)


nk = 1;(*哈密顿量多项式部分展开的阶数*)
Vk = {kx, ky, kz};
varlist = 
  MonomialList[Expand[(1 + Total[Vk])^nk] /. _Integer*b_ -> b, 
   "NegativeDegreeReverseLexicographic"];(*展开多项式并只提取多项式部分，扔掉系数*)
matdim = 4;(*想要展开矩阵的维度*)
numvar = Length[varlist];(*哈密顿量展开,多项式的数量*)
tdc = Array[
   dc, {(Sum[i, {i, matdim}] + Sum[i, {i, matdim - 1}])*numvar}];
matcon1[numvar_, matdim_] := 
 Module[{mdr, mdi}, 
  mdr = Block[{ta, im}, 
    ta = Table[
      IdentityMatrix[matdim], {j, 1, 
       numvar}];(*根据表示的维度构建单位矩阵，数量与多项式的数目相同*)im = 0;
    Do[Do[im = im + 1;
      ta[[i, ii, ij]] = dc[im], {ii, 1, matdim}, {ij, ii, 
       matdim}], {i, 1, numvar}];(*i索引ta保存的每一个单位矩阵，ii和ij则是单位矩阵中元素的索引*)
    Do[Do[ta[[j, ij, ii]] = ta[[j, ii, ij]], {ii, 1, matdim}, {ij, ii,
        matdim}], {j, 1, numvar}];(*这里是对矩阵中的元素做出一些限制,
    也就是令ij位置元素和ji位置元素相同，即对称矩阵要求ij=ji*)ta];
  mdi = Block[{ta, im}, 
    ta = Table[
      IdentityMatrix[matdim], {j, 1, 
       numvar}];(*这里根据多项式展开的数量来确定需要构建的单位矩阵的数量*)
    im = Sum[i, {i, matdim}]*
      numvar;(*这里就是构建一个对称矩阵和反对称矩阵来表示厄米的哈密顿量，而自然这两种矩阵中的参数是要不相同的，所以这里就是用\
im这个索引来区分这些参数，因此在这里的时候使用到的所以就是在上面一个矩阵索引的基础上继续进行叠加，所以这里也就是为什么tdc要构建的dc的\
数量是(Sum[i,{i,id}]+Sum[i,{i,id-1}])*imd个*)
    Do[Do[If[ii == ij, ta[[i, ii, ij]] = 0, im = im + 1;
       ta[[i, ii, ij]] = I*dc[im]], {ii, 1, matdim}, {ij, ii, 
       matdim}], {i, 1, numvar}];(*此时限制对角元素为零，非对角和mdr相同只不过乘以了虚数I*)
    Do[Do[ta[[j, ij, ii]] = -ta[[j, ii, ij]], {ii, 1, matdim}, {ij, 
       ii, matdim}], {j, 1, numvar}];(*此时给出的是限制矩阵为反对称矩阵:ij=-ji*)ta];
  
  {mdr, mdi}](*这里在构建对称矩阵和反对称矩阵的时候使用了一个参量dc来表示矩阵中的参数,用的是全局变量,\
因为后面要对这个参数进行限制*)
matcon[numvar_, matdim_] := 
 Module[{mdr, mdi}, 
  mdr = Block[{ta, im}, 
    ta = Table[
      IdentityMatrix[matdim], {j, 1, 
       numvar}];(*根据表示的维度构建单位矩阵，数量与多项式的数目相同*)im = 0;
    Do[Do[im = im + 1;
      ta[[i, ii, ij]] = tdc[[im]], {ii, 1, matdim}, {ij, ii, 
       matdim}], {i, 1, numvar}];(*i索引ta保存的每一个单位矩阵，ii和ij则是单位矩阵中元素的索引*)
    Do[Do[ta[[j, ij, ii]] = ta[[j, ii, ij]], {ii, 1, matdim}, {ij, ii,
        matdim}], {j, 1, numvar}];(*这里是对矩阵中的元素做出一些限制,
    也就是令ij位置元素和ji位置元素相同，即对称矩阵要求ij=ji*)ta];
  mdi = Block[{ta, im}, 
    ta = Table[
      IdentityMatrix[matdim], {j, 1, 
       numvar}];(*这里根据多项式展开的数量来确定需要构建的单位矩阵的数量*)
    im = Sum[i, {i, matdim}]*
      numvar;(*这里就是构建一个对称矩阵和反对称矩阵来表示厄米的哈密顿量，而自然这两种矩阵中的参数是要不相同的，所以这里就是用\
im这个索引来区分这些参数，因此在这里的时候使用到的所以就是在上面一个矩阵索引的基础上继续进行叠加，所以这里也就是为什么tdc要构建的dc的\
数量是(Sum[i,{i,id}]+Sum[i,{i,id-1}])*imd个*)
    Do[Do[If[ii == ij, ta[[i, ii, ij]] = 0, im = im + 1;
       ta[[i, ii, ij]] = I*tdc[[im]]], {ii, 1, matdim}, {ij, ii, 
       matdim}], {i, 1, numvar}];(*此时限制对角元素为零，非对角和mdr相同只不过乘以了虚数I*)
    Do[Do[ta[[j, ij, ii]] = -ta[[j, ii, ij]], {ii, 1, matdim}, {ij, 
       ii, matdim}], {j, 1, numvar}];(*此时给出的是限制矩阵为反对称矩阵:ij=-ji*)
    ta]; {mdr, mdi}]
hammat[numvar_, matdim_] := 
 Module[{mat1}, mat1 = matcon[numvar, matdim];
  varlist.(mat1[[1]] + mat1[[2]])](*对称矩阵*)
hammatC[numvar_, matdim_] := 
 Module[{mat1}, mat1 = matcon[numvar, matdim];
  varlist.(mat1[[1]] - mat1[[2]])](*反对称矩阵构建*)
diffham[opmat_, opk_, numvar_, matdim_] := 
 opmat.hammat[numvar, matdim].Inverse[
    opmat] - (hammat[numvar, matdim] /. opk)
diffhamC[opmat_, opk_, numvar_, matdim_] := 
 opmat.hammatC[numvar, matdim].Inverse[
    opmat] - (hammat[numvar, matdim] /. opk)
diffcoff[opmat_, opk_, numvar_, matdim_] := 
 Flatten[CoefficientList[diffham[opmat, opk, numvar, matdim], Vk]]
diffcoffC[opmat_, opk_, numvar_, matdim_] := 
 Flatten[CoefficientList[diffhamC[opmat, opk, numvar, matdim], Vk]]
constrain[opmat_, opk_, numvar_, matdim_] := 
 Module[{caplist(*这个参数我是想作为方程求解的时候的变量,因为通过对称性限制之后,一些参数可能会被限制为零,
   那么在方程求解过程中这些参量就不会作为变量来求解了*), opcons(*这个参量存储对称限制后相差的部分*), 
   re1(*这个变量就用来存储系数方程的解*)},
  opcons = diffcoff[opmat, opk, numvar, matdim];
  re1 = Flatten[
    Quiet@Solve[opcons == Table[0, {i, 1, Length[opcons]}]]]; 
  tdc = (tdc /. re1)]


constrainC[opmat_, opk_, numvar_, matdim_] := 
 Module[{caplist(*这个参数我是想作为方程求解的时候的变量,因为通过对称性限制之后,一些参数可能会被限制为零,
   那么在方程求解过程中这些参量就不会作为变量来求解了*), opcons(*这个参量存储对称限制后相差的部分*), 
   re1(*这个变量就用来存储系数方程的解*)},
  opcons = diffcoffC[opmat, opk, numvar, matdim];
  re1 = Flatten[
    Quiet@Solve[opcons == Table[0, {i, 1, Length[opcons]}]]]; 
  tdc = (tdc /. re1)]

constrainU[OpRes_, numvar_, matdim_] := 
 Module[{coff}, 
  Do[coff = 
    constrain[OpRes[[i, 1]], OpRes[[i, 2]], numvar, matdim], {i, 1, 
    Length[OpRes]}];
  coff](*遍历所有的幺正操作操作进行约束*)


constrainAU[OpRes_, numvar_, matdim_] := 
 Module[{coff}, 
  Do[coff = 
    constrainC[OpRes[[i, 1]], OpRes[[i, 2]], numvar, matdim], {i, 1, 
    Length[OpRes]}];
  coff](*遍历所有的反幺正操作操作进行约束*)


constrainAll[OpRes_, numvar_, matdim_] := 
 Module[{}, constrainU[OpRes[[1]], numvar, matdim]; 
  constrainAU[OpRes[[2]], numvar, matdim];](*整合幺正与非幺正操作*)


constrainAll[OperationRes, numvar, matdim]
tdc = Flatten[
   tdc /. Quiet[
     Solve[Variables[tdc] == 
       Table[Subscript[c, i], {i, 1, Length[Variables[tdc]]}], 
      Variables[tdc]]]];(*这里看起来是个方程求解,其实就是想将之前的构建矩阵时候的参量换个形式表示出来*)
Style[MatrixForm[
  FullSimplify[
   hammat[numvar, matdim] /. {kx -> Subscript[k, x], 
     ky -> Subscript[k, y], kz -> Subscript[k, z]}]], 
 FontFamily -> "Times", 20, Blue]
```