#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import signal as sgn
import matplotlib.pyplot as plt
import control

def fDegrau(t_array,tamanho=1):
    y=[]
    for t in t_array:
        tmp=0
        if t>=0:
            tmp=tamanho
        y.append(tmp)
    return y

def getLastStableEl(array):
    lastValue=None
    count=0
    last_count=0
    last_lastValue=None
    for el in array:
        if not lastValue:
            lastValue=el
        else:
            if el!=lastValue:
                last_lastValue=lastValue
                last_count=count
                count=0
                lastValue=el
            else:
                count+=1

    if count>=last_count:
        out=lastValue
    else:
        out=last_lastValue
    return out
   
def getGain(array,input_value):
    return getLastStableEl(array)/input_value

def getTimeConst(t_array,array):
    last=getLastStableEl(array)
    for i in range(len(t_array)):
        if array[i]>last*0.632:
            return t_array[i]
        

# a
# calcular

# b
t_inicial=0
t_final=10
passo_T=0.001
time_array = np.arange(t_inicial, t_final,passo_T)
tfs = control.TransferFunction([.5],[.5,1])
print('Ganho DC real:',tfs.dcgain())
mecum = control.forced_response(tfs,time_array,fDegrau(time_array))

print('Ganho DC:',getGain(mecum[1],input_value=1))
print('Constante de tempo:',getTimeConst(mecum[0],mecum[1]))

plt.plot(mecum[0],mecum[1],label='Resposta ao degrau')
plt.legend(loc='best')
plt.show()


# c
T=0.1
time_array_disc = np.arange(t_inicial, t_final,T)
tfz = control.c2d(tfs, T)
print(tfz)

#d
print('Ganho DC real:',tfz.dcgain())
cum = control.forced_response(tfz,time_array_disc,fDegrau(time_array_disc))
print('Ganho DC:',getGain(cum[1],input_value=1))
print('Constante de tempo:',getTimeConst(cum[0],cum[1]))

#e
plt.step(cum[0],cum[1],label='Resposta ao degrau discreto')
plt.plot(mecum[0],mecum[1],label='Resposta ao degrau continuo')
plt.legend(loc='best')
plt.show()


#f
k_p=2
tfs_MF = (k_p*tfs)/(1+k_p*tfs)
print(tfs_MF)

#g
# calcular


#h
print('Ganho DC real:',tfs_MF.dcgain())
feden = control.forced_response(tfs_MF,time_array,fDegrau(time_array,tamanho=100))
print('Ganho DC:',getGain(feden[1],input_value=100))
print('Constante de tempo:',getTimeConst(feden[0],feden[1]))
plt.plot(feden[0],feden[1],label='Resposta ao degrau continuo MF')
plt.legend(loc='best')
plt.show()


#i 
k_p=2
tfz_MF = (k_p*tfz)/(1+k_p*tfz)
print(tfz_MF)

#j
# calcular

#k
print('Ganho DC real:',tfz_MF.dcgain())
xunda = control.forced_response(tfz_MF,time_array_disc,fDegrau(time_array_disc,tamanho=100))
print('Ganho DC:',getGain(xunda[1],input_value=100))
print('Constante de tempo:',getTimeConst(xunda[0],xunda[1]))
plt.plot(feden[0],feden[1],label='Resposta ao degrau continuo MF')
plt.step(xunda[0],xunda[1],label='Resposta ao degrau discreto MF')
plt.legend(loc='best')
plt.show()

#l
kps=[2,4,6,8,10,12,14,16,18,20]
for k_p in kps:
    tfs_MF = (k_p*tfs)/(1+k_p*tfs)
    bund = control.forced_response(tfs_MF,time_array,fDegrau(time_array,tamanho=100))
    plt.plot(bund[0],bund[1],label='Resposta ao degrau continuo MF, Kp={}'.format(k_p))
plt.legend(loc='best')
plt.show()

#m 
kps.reverse()
for k_p in kps:
    tfz_MF = (k_p*tfz)/(1+k_p*tfz)
    zeleoterio = control.forced_response(tfz_MF,time_array_disc,fDegrau(time_array_disc,tamanho=100))
    plt.step(zeleoterio[0],zeleoterio[1],label='Resposta ao degrau discreto MF, Kp={}'.format(k_p))
plt.legend(loc='best')
plt.show()

