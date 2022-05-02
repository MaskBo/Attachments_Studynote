# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:06:20 2022

@author: 14538
"""

import math
import matplotlib.pyplot as plt

#定义常数
N = 100 #积分精度
pi = 3.1415 #常数pi
me = 0.51 #电子质量，单位1MeV
MX = [0.1,0.5,1.0] #DPs的质量，单位MeV，后面会更具画图重新改变
alpha = 197/137 #精细结构常数,单位 10^-13  MeV*cm
gx2 = 4*pi*197/137 #暗光子相互作用顶角,单位 10^-13  MeV*cm
Zls = 3 #液体闪烁体平均原子量

Ne = 3.51 #探测器内液体闪烁体LS的密度10^29 每吨
R = 3 #探测器距离反应堆核心的距离，单位10^3cm
T = 1 #探测器曝光时间


#类逆康普顿散射过程的微分截面
def DCDEr(Ex,Er):
    D1 = 16*Ex**4*me**3-16*Ex**2*me**3*Er*(2*Ex+me)+2*me*mx**4*(3*Ex**2-5*Ex*Er+me**2+2*Er*(Er+me))
    D2 = 8*me**3*Er**2*(3*Ex**2+2*Ex*me+me**2)+4*me**2*mx**2*(4*Ex**3-8*Ex**2*Er+Ex*Er*(5*Er-2*me))
    D3 = 4*me**2*mx**2*Er*(2*me**2+3*me*Er-Er**2)+mx**6*(Ex+me-Er)-8*Ex*me*3*Er**3
    D = D1+D2+D3
#    print(D)
    dcder = Zls*alpha*gx2/(8*me**2*(Ex**2-mx**2)*(Ex-Er)**2*(2*Ex*me+mx**2)**2)*(D)
#    print(dcder)
    return dcder;

#对于特定能量Ex的暗物质生成光子能量范围
def ErM(Ex):
    ErMin = (0.5*mx**2+Ex*me)/(me+Ex+math.sqrt(Ex**2-mx**2))
    ErMax = (0.5*mx**2+Ex*me)/(me+Ex-math.sqrt(Ex**2-mx**2))
    return ErMin,ErMax;

#计算特定能量暗物质Ex的总散射
def ToTalCS(Ex):
    Emin,Emax = ErM(Ex)
    Deta = (Emax-Emin)/N
#    print(Ex,Emin,Emax)
    totalcs = 0
    for i in range(N):
        Er = Emin + i*Deta
        totalcs = totalcs + DCDEr(Ex,Er)*Deta
    print(Ex,totalcs)
    return totalcs;

plt.figure(figsize=(6.5,5),dpi=330)
Line_Name = []
Line_Flag = 0

for mx in MX:
    filename = './datas/DataEx_DN2_'+str(mx)+'.txt'
    # 读取反应堆产生暗物质流的数据
    infile = open(filename,'r')
    
    Ex = [] #需要探测的暗物质能量
    DN2DE2 = [] #暗物质流密度
    for line in infile:
        ab = line.split()
        Ex.append(float(ab[0]))
        DN2DE2.append(float(ab[1]))
    infile.close()
    
    #探测效率
    Exnew = []
    DNobs = []
    for i in range(len(Ex)):
        ex = Ex[i]+0.001
#        if Ex[i]<0.8:
#            continue
        Exnew.append(ex)
        DNobs.append(-ToTalCS(ex)*DN2DE2[i]*Ne/(4*pi*R**2))
    Line_Name.append('$m_X =$'+str(mx)+'MeV')
    if Line_Flag == 0:
        Line1, = plt.plot(Exnew,DNobs)
    if Line_Flag == 1:
        Line2, = plt.plot(Exnew,DNobs)
    if Line_Flag == 2:
        Line3, = plt.plot(Exnew,DNobs)
    Line_Flag +=1



plt.yscale('log')
plt.xlabel('$E_X$(MeV)')
plt.ylabel('$\\frac{d N_{obs}}{d E_X},\\frac{10^{18}}{MeV\cdot s}$')
plt.xlim(0,5)
plt.ylim(10**(-2),10**4)
plt.legend((Line1,Line2,Line3),Line_Name)
plt.savefig('TAO_Detection_Rate.png')
plt.show()
