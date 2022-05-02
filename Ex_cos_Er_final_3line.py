# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 10:20:10 2022

@author: 14538
"""

import math
import matplotlib.pyplot as plt

me = 0.51
Mx = [0.1,0.5,1.0,1.5]


def equation1():
    result = (-0.5*mx**2+Ex*me)/(me-Ex+math.sqrt(Ex**2-mx**2)*co)
    return result;


M = 50
N = 3000

EEx = [0]*M


for mx in Mx:


    fig = plt.figure(figsize=(7,5.5),dpi=330)
    ax = plt.axes(projection='3d')
    
    for i in range(M):
        Ex = 0.05*i+mx         #入射光子能量MeV
        EEx[i] = Ex
        theta = []
        er1 = []
        for j in range(N):
            co = 2*j/(N-1)-1    #cos(theta)的取值
    
            er_possible = equation1()
            if er_possible<0:
                continue
            if er_possible<Ex-me:
                continue
            if er_possible>200:
                er_possible = 200
            er1.append(er_possible)
            theta.append(co)
    
        x = theta
        y = EEx
        z = er1
        if i==0:
            continue
        if x == [] or z == []:
            continue
        ax.plot3D([y[i]]*len(x),x,z,'gray')  
    
    #ax.set_zscale('log')
    ax.set_ylabel('cos$\Theta$')
    ax.set_zlabel('$E_\gamma$(MeV)')
    ax.set_xlabel('$E_X$(MeV)')
    ax.set_xlim(mx,mx+0.05*M)
    ax.set_zlim(0,200)
    
    filename = './figs/Ex_cos_Er_'+str(mx)+'.png'
    plt.savefig(filename)
    plt.show