# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 23:03:56 2020

@author: zhangshuo
"""


from qutip import *
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import time
 
 
start = time.time()


#Hilbert dimension
N = 20
# parameters
nu = 1
nbar = 3#初始声子数
gammag = 40/3
gammar = 20/3
gamma = gammag + gammar
Omegag = 10
Omegar = 20
deltag = (Omegag*Omegag+Omegar*Omegar)/(4*nu)-nu#最优条件下的失谐量
deltar = (Omegag*Omegag+Omegar*Omegar)/(4*nu)-nu
deltagr = deltag - deltar
ldg = 0.15
ldr = -0.15
ld = ldg - ldr
Omegaeff = ld*Omegag*Omegar/(2*deltag)

#T1和T2
taus1 = np.linspace(0, 5, 11)
taus2 = np.linspace(0, 0.3, 2)

maxt = taus1[-1]+taus2[-1]#T1+T2时长
maxt2 = taus2[-1]#T2时长


#operators
gx1 =Qobj([[0, 0, 1],[0, 0, 0],[1, 0, 0]])
gy1 =Qobj([[0, 0, -1j],[0, 0, 0],[1j, 0, 0]])
rx1 =Qobj([[0, 0, 0],[0, 0, 1],[0, 1, 0]])
ry1 =Qobj([[0, 0, 0],[0, 0, -1j],[0, 1j, 0]])
grx1 =Qobj([[0, 1, 0],[1, 0, 0],[0, 0, 0]])
gry1 =Qobj([[0, -1j, 0],[1j, 0, 0],[0, 0, 0]])
ee1 =Qobj([[0, 0, 0],[0, 0, 0],[0, 0, 1]])
rr1 =Qobj([[0, 0, 0],[0, 1, 0],[0, 0, 0]])
gg1 =Qobj([[1, 0, 0],[0, 0, 0],[0, 0, 0]])
gm1 =Qobj([[0, 0, 1],[0, 0, 0],[0, 0, 0]])
rm1 =Qobj([[0, 0, 0],[0, 0, 1],[0, 0, 0]])
grm1 =Qobj([[0, 1, 0],[0, 0, 0],[0, 0, 0]])
atom0 = Qobj([[0, 0, 0],[0, 0, 0],[0, 0, 0]])

a = tensor(destroy(N), qeye(3))
gx = tensor(qeye(N), gx1)
gy = tensor(qeye(N), gy1)
rx = tensor(qeye(N), rx1)
ry = tensor(qeye(N), ry1)
grx = tensor(qeye(N), grx1)
gry = tensor(qeye(N), gry1)
ee = tensor(qeye(N), ee1)
rr = tensor(qeye(N), rr1)
gg = tensor(qeye(N), gg1)
gm = tensor(qeye(N), gm1)
rm = tensor(qeye(N), rm1)
grm = tensor(qeye(N), grm1)
idd = tensor(qeye(N), qeye(3))

#time-independent Hamiltonian
H00= -1*deltag * ee+nu* (a.dag()* a)-( 1* deltagr * rr)
H01= 0.5*Omegag*(gm*(-1j*ldg*(a+a.dag())).expm()+gm.dag()*(1j*ldg*(a+a.dag())).expm())
H02= 0.5*Omegar*(rm*(-1j*ldr*(a+a.dag())).expm()+rm.dag()*(1j*ldr*(a+a.dag())).expm())

#第一段演化的哈密顿
H1 = H00 + H01 + H02

#第二段演化的哈密顿 
H2 = nu* (a.dag()* a)+ H02

# 耗散
c_ops=[]
mmm = np.linspace(-1, 1, 129)
for i in mmm:
    xxx =np.sqrt(gammag*(1+ i*i)/172.6719) * gm*(-1j*ldg*(a+a.dag())*i).expm()
    yyy =np.sqrt(gammar*(1+ i*i)/172.6719) * rm*(-1j*ldr*(a+a.dag())*i).expm()
    c_ops.append(xxx)
    c_ops.append(yyy)




#初始声子态为热态，内态为暗态
darkstate = 1/np.sqrt(Omegag*Omegag+Omegar*Omegar)*(Omegar*basis(3,0)-Omegag*basis(3,1))
rho00 = darkstate*darkstate.dag()
rho00 = tensor(thermal_dm(N,nbar),rho00)#初态



#演化开始
rho0 = rho00
rhos = [1] #用以记录不同时刻的密度算子
taus = np.array([0]) #用以记录总时间
#result
aaaa = range(10) #周期数
for i in aaaa: #一个周期的演化
    # 第一段演化
    print(i+1)
    del rhos[-1] #本段初态与上一段末态重复，需删除一个
    data1 = mesolve(H1, rho0, taus1, c_ops, [ ]) #T1的演化
    taus = np.append(taus,maxt*(i)+taus1) #T1加入taus
    taus = np.delete(taus, -1*len(taus1)) #本段末时刻与下一段初时刻重复，需删除一个
    rhos = rhos + data1.states  #T1的演化的态加入rhos
    #第二段演化
    rho0 = rhos[-1] #第二段的初态
    del rhos[-1] #本段初态与上一段末态重复，需删除一个
    data2= mesolve(H2, rho0, taus2, c_ops, [ ]) #T2的演化
    rhos = rhos + data2.states #T2的演化的态加入rhos
    taus = np.append(taus,(maxt*(i+1)-maxt2)+taus2) #T2加入taus
    taus = np.delete(taus, -1*len(taus2)) #本段末时刻与下一段初时刻重复，需删除一个
    rho0 = rhos[-1] #下一段的初态   
    
#如需继续演化激活这一段
'''
print('开始继续演化')
del rhos[-1] #删除重复态
taus = np.delete(taus, -1) #删除重复时间
taus3 = linspace(0, 500, 1000)
taus = np.append(taus,maxt*(i+1)+taus3)
data3 = mesolve(H1, rho0, taus3, c_ops, [ ])
rhos = rhos + data3.states
'''

#taus时间段的对比演化
print('开始对比演化')
data4 = mesolve(H1, rho00, taus, c_ops, [ ])


#理论参数
gammad = (gammag*(Omegar*Omegar)+gammar*(Omegag*Omegag))/(Omegag*Omegag+Omegar*Omegar)
gammageff=(gammad)*nu/sqrt(deltag*deltag+Omegag*Omegag+Omegar*Omegar)#等效耗散
www = gammageff/2*(1/(nbar+1))#强耦合冷却速度的上限
n1=nbar-1/2*(nbar/(nbar+1))#强耦合下的初始声子数
Omegaeff = Omegag*Omegar/sqrt(Omegag*Omegag+Omegar*Omegar) #omegad
nss = gamma*gamma/(16*deltag*deltag)+(ld*ld*Omegaeff*Omegaeff)*(1/deltag)*(1/8) #理论的EIT最终声子数

plt.figure(1)#创建图表1  
nnn1 = expect(a.dag()*a, rhos)         #经调控的声子数
nnn2 = expect(a.dag()*a, data4.states) #原EIT声子数
plt.semilogy(taus,nnn2,color="blue",linewidth=2.5) #原EIT冷却曲线，蓝色
plt.semilogy(taus,nnn1,color="red",linewidth=2.5)  #经调控的冷却曲线，红色
plt.semilogy(taus,n1*e**(-www*taus)+nss,color="black",linewidth=2.5)  #强耦合冷却速度上限，黑色
ax=plt.gca()
ax.set_ylabel('<n>',fontsize=22,labelpad = 1)
ax.set_xlabel('t',fontsize=22,labelpad = 0)
ax.spines['bottom'].set_linewidth(1.7)

ax.spines['top'].set_linewidth(1.7)

ax.spines['left'].set_linewidth(1.7)

ax.spines['right'].set_linewidth(1.7)
end = time.time()

print("Execution Time: ", end - start)
show()