# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 10:48:07 2022

@author: 14538
"""

import math
import matplotlib.pyplot as plt

#===================================定义常数=================================


alpha = 1/137*19.7;          #精细结构常数
Z = 90;                 #原子序数
g = 1/137*19.7*4*3.14;              #DPs的耦合系数

me = 0.51;              #单位MeV
MX = [0.1,0.5,1.0];     #DPs的质量单位MeV
P = 1;                  #反应堆功率GW



#===================================定义函数=================================


def s(E1):     
    '''
     Mandelsm变量s的定义
    '''
    s = me**2+2*E1*me 
    return s;

def E1M(E2):
    '''
    要打出E2能量的DPs需要入射光子能量下限
    '''
    e1m =(-0.5*mx**2+E2*me)/(me-E2 +math.sqrt(E2**2-mx**2))
    return e1m;

def x(E1,E2):
    '''
    分数动量定义
    '''
    x = 1- E2/E1 + mx**2/(2*E1*me)
    return x;

def DN1dE1(E1):
    '''
    由反应堆特性决定的光子流密度，单位10^21
    '''
    DN = 0.58*P*math.exp(-1.1*E1)
    return DN;

def DifCS(E1,x,s):
    '''
    计算微分散射截面
    '''
    C1 = -2*me**2*s*(x**2+2) + 2*mx**4 - 2*mx**2*s*x + (x**2-2*x+2)*s**2 \
    -(me**4*(x**3-3*x**2-2)+2*me**2*mx**2*(x-2))/(1-x)
    DCS = 1/E1*Z*alpha*g**2/(2*(s-me**2)**3*(1-x))*C1
#    DCS = 1/E1*Z*alpha*g**2/(2*(s-me**2)**3*(1-x))*C1
    return DCS;

def TotalCS(s):
    '''
    计算总反应截面
    '''
    A1 = me**6 -me**4*(s+mx**2) +me**2*s*(2*mx**2+15*s) +s**2*(s+7*mx**2)
    A = math.sqrt(me**4-2*me**2*(s+mx**2)+(s-mx**2)**2)*A1
    B = -2*s**2*(-3*me**4+2*me**2*(mx**2-3*s)+2*mx**4-2*mx**2*s+s**2)
    c1 = me**2-math.sqrt(me**4-2*me**2*(mx**2+s)+(s-mx**2)**2)-mx**2+s
    c2 = me**2+math.sqrt(me**4-2*me**2*(mx**2+s)+(s-mx**2)**2)-mx**2+s
    C = c1/c2
#    print(A,B,C)
#    TCS = Z*alpha*g**2/(8*s**2*(s-me**2)**3)*(A+B*math.log(C))*math.sqrt(e2**2-mx**2)
    TCS = Z*alpha*g**2/(4*s**2*(s-me**2)**3)*(A+B*math.log(C))
    return TCS;



#===============由于光子与核反应堆的散射截面由列表给出下面读取===================


infile = open('./datas/totaldata.txt','r')         #打开文件

E_gamma_data = []                                           #存放光子能量
PhotonTTCS = []                                    #光子与反应堆物质反应的总截面
for line in infile:                                #按行读取
    ab = line.split()                              #每行数据按空格分开
    E_gamma_data.append(float(ab[0]))                       #读取数据
    PhotonTTCS.append(float(ab[1]))
infile.close()                                     #关闭文件



#============================主函数部分======================================

plt.figure(figsize=(6.8,5),dpi=330)
Line_Name = []
line_flag = 0

for mx in MX:
    E2 = []          #定义出射暗物质能量,作图的横坐标
    DN2DE2 = []      #定义光子产生率
    N = 200          #将黎曼积分分成若干份
    M = 500          #将Ex分成M份，即画M个点
    
    if mx>me:  #产生暗物质的最低能量
        E2Min = me/2+mx**2/(2*me)+0.00001
    else:
        E2Min = mx
             
    for i in range(M):
        e2 = E2Min+i*(5.0-E2Min)/M
        E2.append(e2)
        E1Min = E1M(e2)
        E1Max = 20
        dN2dE2 = 0
        for j in range(N):
            Deta = (E1Max-E1Min)/N
            e1 = E1Min+j*Deta
            k = 0
            for e11 in E_gamma_data:   #从列表中抽出对应光子能量e1的中反应截面
#                if e11<0.8 or e11>10:
#                    PhotonTotalCS = 0
                if e1<e11:
                    break
                k = k+1
            PhotonTotalCS = PhotonTTCS[k-1]
#            PhotonTotalCS = PhotonTTCS[k-1]/(19.7**2*4*3.14)
            dN2dE2 = dN2dE2 + 1/(PhotonTotalCS+TotalCS(s(e1)))*DifCS(e1,x(e1,e2),s(e1))*DN1dE1(e1)*Deta
#            if e2<1.4 and j==0:
#                print(e1,e2,PhotonTotalCS,TotalCS(s(e1)),DifCS(e1,x(e1,e2),s(e1))*Deta,DN1dE1(e1))
        DN2DE2.append(dN2dE2)
    Line_Name.append('$m_X =$'+str(mx)+'MeV')
    if line_flag == 0:
            Line1, = plt.plot(E2,DN2DE2)
    if line_flag == 1:
            Line2, = plt.plot(E2,DN2DE2)
    if line_flag == 2:
            Line3, = plt.plot(E2,DN2DE2)
    line_flag +=1




#=========================将计算数据写入文件===============================

    
    filename = './datas/DataEX_DN2_' + str(mx)+'.txt'    #定义文件名
    outfile = open(filename, 'w')                #以写入打开文件
    for i in range(len(E2)):                     #循环写入数据
        outfile.write('%14.8f %14.8f\n' %(E2[i],DN2DE2[i]))
    outfile.close()                              #写入完成并关闭文件



#===============================绘图参数==================================

plt.yscale('log')
plt.ylim(10**(-3),10**(1))
plt.xlabel('$E_X$(MeV)')
plt.ylabel('$\\frac{dN_X}{dE_X},\\frac{10^{21}}{MeV\cdot s}$')
plt.legend((Line1,Line2,Line3),Line_Name)
plt.savefig('TAO_Product_Rate.png')
plt.show()
