---
title: 超大型系数矩阵对角化
tags: Python Study
layout: article
license: true
key: a20200810
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
toc: true
pageview: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---


最近在进行一些三维实空间中系统性质的研究，利用紧束缚束缚模型，当三个方向都取开边界时，矩阵会变的非常大，组里的服务器完全不能计算，所以只好寻找一些替代的方法。首先哈密顿量是个大型系数矩阵，这个时候可以通过[lanczos 算法](https://en.wikipedia.org/wiki/Lanczos_algorithm)将对称的厄密矩阵变成三对角矩阵，这个时候矩阵的本征值求解就变的比较容易了。
<!--more-->
# 算法介绍
```python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
# import useful modules
import matplotlib 
from math import factorial
from numpy import *
from pylab import *
from numpy.polynomial.hermite import *
 
# use LaTeX, choose nice some looking fonts and tweak some settings
matplotlib.rc('font', family='serif')
matplotlib.rc('font', size=16)
matplotlib.rc('legend', fontsize=16)
matplotlib.rc('legend', numpoints=1)
matplotlib.rc('legend', handlelength=1.5)
matplotlib.rc('legend', frameon=False)
matplotlib.rc('xtick.major', pad=7)
matplotlib.rc('xtick.minor', pad=7)
matplotlib.rc('lines', lw=1.5)
matplotlib.rc('text', usetex=True)
matplotlib.rc('text.latex', 
              preamble=[r'\usepackage[T1]{fontenc}',
                        r'\usepackage{amsmath}',
                        r'\usepackage{txfonts}',
                        r'\usepackage{textcomp}'])
 
close('all')
figure(figsize=(6, 4.5))
 
c=137.0359998  # speed of light in a.u.
N=1024+512
 
# the 1d Dirac Hamiltonian
def H(Psi, x, t, dt):
    Psi=reshape(Psi, (N, 2))
    dx=x[1]-x[0]
    Psi_new=empty_like(Psi)
    Psi_new[1:-1, 0]=-1j*c*(Psi[2:, 1]-Psi[0:-2, 1])/(2*dx) + c**2*Psi[1:-1, 0]
    Psi_new[1:-1, 1]=-1j*c*(Psi[2:, 0]-Psi[0:-2, 0])/(2*dx) - c**2*Psi[1:-1, 1]
    Psi_new[0, 0]=Psi_new[0, 1]=0
    Psi_new[-1, 0]=Psi_new[-1, 1]=0
    Psi_new*=dt
    return reshape(Psi_new, 2*N)
 
# the Lanczos algorithm
def Lanczos(Psi, x, t, dt, H, m):
    Psi_=Psi.copy()
    # run Lanczos algorithm to calculate basis of Krylov space
    V_j, V_j_1=[], []
    A=zeros((m, m))
    norms=zeros(m)
    for j in range(0, m):
        norms[j]=norm(Psi_)
        V_j=Psi_/norms[j]
        Psi_=H(V_j, x, t, dt)
        if j>0:
            A[j-1, j]=A[j, j-1]=vdot(V_j_1, Psi_).real
            Psi_-=A[j-1, j]*V_j_1
        A[j, j]=vdot(V_j, Psi_).real
        Psi_-=A[j, j]*V_j
        V_j_1, V_j=V_j, V_j_1
    # diagonalize A
    l, v=eig(A)
    # calculate matrix exponential in Krylov space
    c=dot(v, dot(diag(exp(-1j*l)), v[0, :]*norms[0]))
    # run Lanczos algorithm 2nd time to transform result into original space
    Psi_, Psi=Psi, zeros_like(Psi)
    for j in range(0, m):
        V_j=Psi_/norms[j]
        Psi+=c[j]*V_j
        Psi_=H(V_j, x, t, dt)
        if j>0:
            A[j-1, j]=A[j, j-1]=vdot(V_j_1, Psi_).real
            Psi_-=A[j-1, j]*V_j_1
        A[j, j]=vdot(V_j, Psi_).real
        Psi_-=A[j, j]*V_j
        V_j_1, V_j=V_j, V_j_1
    return Psi
 
# define computational grid
x0, x1=-0.5, 0.5
x=linspace(x0, x1, N)
dx=x[1]-x[0]       # size of spatial grid spacing
dt=4./c**2         # temporal step size
 
# constuct momentum grid
dp=2*pi/(N*dx)
p=(arange(0, N)-0.5*(N-1))*dp
# choose initial condition
p_mu=75.       # mean momentum
sigma_p=50.    # momentum width
x_mu=-0.05     # mean position
# upper and lower components of free particle states, 
# see e.g. Thaller »Advanced visual quantum mechanics«
d_p=sqrt(0.5+0.5/sqrt(1+p**2/c**2))
d_m=sqrt(0.5-0.5/sqrt(1+p**2/c**2))
d_m[p<0]*=-1
# initial condition in momentum space, 
# gaussian wave packet of positive energy states
rho=(2*pi*sigma_p**2)**(-0.25)*exp(-(p-p_mu)**2/(4*sigma_p**2) - 1j*p*x_mu) 
Psi=zeros((N, 2), dtype='complex')
Psi[:, 0]=d_p*rho
Psi[:, 1]=d_m*rho
# transform into real space with correct complex phases 
Psi[:, 0]*=exp(1j*x[0]*dp*arange(0, N))
Psi[:, 1]*=exp(1j*x[0]*dp*arange(0, N))
Psi=ifft(Psi, axis=0)
Psi[:, 0]*=dp*N/sqrt(2*pi)*exp(1j*x*p[0])
Psi[:, 1]*=dp*N/sqrt(2*pi)*exp(1j*x*p[0])
 
# propagate 
for k in range(0, 20):
    # plot wave function
    clf()
    plot(x, Psi[:, 0].real**2+Psi[:, 0].imag**2+
         Psi[:, 1].real**2+Psi[:, 1].imag**2, 
         color='#266bbd', label=r'$|\Psi(x)|^2$')
    gca().set_xlim(x0, x1)
    gca().set_ylim(-1, 16)
    xlabel(r'$x$')
    legend(loc='upper left')
    tight_layout()
    draw()
    show()
    pause(0.05)
    Psi=reshape(Lanczos(reshape(Psi, 2*N), x, 0, dt, H, 128), (N, 2))
```

# 代码展示
```python
import numpy as np
from spartan import expr, util
import math

def solve(A, AT, desired_rank, is_symmetric=False):
  '''
  A simple implementation of the Lanczos algorithm
  (http://en.wikipedia.org/wiki/Lanczos_algorithm) for eigenvalue computation.
  Like the Mahout implementation, only the matrix*vector step is parallelized.
  
  First we use lanczos method to turn the matrix into tridiagonoal form. Then
  we use numpy.linalg.eig function to extract the eigenvalues and eigenvectors 
  from the tridiagnonal matrix(desired_rank*desired_rank). Since desired_rank 
  should be smaller than the size of matrix, so we could it in local machine 
  efficiently. 
  '''
  # Calculate two more eigenvalues, but we only keep the largest desired_rank
  # one. Doing this to keep the result consistent with scipy.sparse.linalg.svds.
  desired_rank += 2

  n = A.shape[1]
  v_next = np.ones(n) / np.sqrt(n)
  v_prev = np.zeros(n)
  beta = np.zeros(desired_rank+1)
  beta[0] = 0
  alpha = np.zeros(desired_rank)

  # Since the disiredRank << size of matrix, so we keep
  # V in local memory for efficiency reason(It needs to be updated
  # for every iteration). 
  # If the case which V can't be fit in local memory occurs, 
  # you could turn it into spartan distributed array. 
  V = np.zeros((n, desired_rank))


  for i in range(0, desired_rank):
    util.log_info("Iter : %s", i)
    v_next_expr = expr.from_numpy(v_next.reshape(n, 1))

    if is_symmetric:
      w = expr.dot(A, v_next_expr).optimized().glom().reshape(n)
    else:
      w = expr.dot(A, v_next_expr)
      w = expr.dot(AT, w).optimized().glom().reshape(n)

    alpha[i] = np.dot(w, v_next)
    w = w - alpha[i] * v_next - beta[i] * v_prev
    
    # Orthogonalize:
    for t in range(i):
      tmpa = np.dot(w, V[:, t])
      if tmpa == 0.0:
        continue
      w -= tmpa * V[:, t] 

    beta[i+1] = np.linalg.norm(w, 2) 
    v_prev = v_next
    v_next = w / beta[i+1]
    V[:, i] = v_prev
  
  # Create tridiag matrix with size (desired_rank X desired_rank)  
  tridiag = np.diag(alpha)
  for i in range(0, desired_rank-1):
    tridiag[i, i+1] = beta[i+1] 
    tridiag[i+1, i] = beta[i+1]
  
  # Get eigenvectors and eigenvalues of this tridiagonal matrix.  
  # The eigenvalues of this tridiagnoal matrix equals to the eigenvalues
  # of matrix dot(A, A.T.). We can get the eigenvectors of dot(A, A.T) 
  # by multiplying V with eigenvectors of this tridiagonal matrix.
  d, v = np.linalg.eig(tridiag) 
  
  # Sort eigenvalues and their corresponding eigenvectors 
  sorted_idx = np.argsort(np.absolute(d))[::-1]
  d = d[sorted_idx]
  v = v[:, sorted_idx]
  
  # Get the eigenvetors of dot(A, A.T)
  s = np.dot(V, v)
  return d[0:desired_rank-2], s[:, 0:desired_rank-2] 
```
# Ref
[code1](https://web.cs.ucdavis.edu/~bai/ET/lanczos_methods/lanczos_methods.html)
