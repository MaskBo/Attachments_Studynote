# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 20:44:35 2022

@author: 14538
"""

import math
import matplotlib.pyplot as plt

me = 0.51
mx = 0.1


def Sx(Er):      #分数动量定义
    s = me**2+2*Er*me
    return s;

def equation1(s):
    result = math.sqrt(co**2*(s-me**2)**2*(me**6-s**2*si**2*mx**2+me**4*(-2*s+(-3+co**2)*mx**2)+me**2*(s**2-2*s*co**2*mx**2+mx**4)))
    return result;

def equation2(s):
    eq2 = s**2*me-me**5+s*me*mx**2+me**3*mx**2
    return eq2;

def equation3(s):
    eq3 = s**2*si**2+s*(3+co**2-si**2)*me**2+si**2*me**4
    return eq3;

M = 10
N = 1000
theta = [0]*N
EEr = [0]*M
ex1 = [0]*N
ex2 = [0]*N


fig = plt.figure(figsize=(8,6.5),dpi=330)
ax = plt.axes(projection='3d')

for i in range(M):
    Er = 0.5*i+0.5         #入射光子能量MeV
    EEr[i] = Er
    for j in range(N):
        s = Sx(Er)
        co = 2*j/(N-1)-1    #cos(theta)的取值
        si = math.sqrt(1-co**2)
        theta[j] = co
        ex1[j]=(equation2(s)-equation1(s))/equation3(s)
        if j == 99:
            print((equation2(s)+equation1(s))/equation3(s))
        ex2[j]=(equation2(s)+equation1(s))/equation3(s)
    x = theta
    y = EEr
    z1 = ex1
    z2 = ex2
    ax.plot3D(x,[y[i]]*N,z1,'gray')
#    ax.plot3D(x,[y[i]]*N,z2,'green')

for i in range(M):
    Er = 0.5*i+0.5         #入射光子能量MeV
    EEr[i] = Er
    for j in range(N):
        s = Sx(Er)
        co = 2*j/(N-1)-1    #cos(theta)的取值
        si = math.sqrt(1-co**2)
        theta[j] = co
        ex1[j]=(equation2(s)-equation1(s))/equation3(s)
        if j == 99:
            print((equation2(s)+equation1(s))/equation3(s))
        ex2[j]=(equation2(s)+equation1(s))/equation3(s)
    x = theta
    y = EEr
    z1 = ex1
    z2 = ex2
#    ax.plot3D(x,[y[i]]*N,z1,'gray')
    ax.plot3D(x,[y[i]]*N,z2,'green')


ax.set_xlabel('cos$\Theta$')
ax.set_ylabel('$E_\gamma$(MeV)')
ax.set_zlabel('$E_X$(MeV)')

plt.savefig('./figs/Cos_Ex_Er.png')
plt.show
