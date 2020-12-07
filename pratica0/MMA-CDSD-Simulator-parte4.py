#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import signal as sgn
import matplotlib.pyplot as plt
import control

def closestMin(array,min_ref):
    array=array.copy()
    array.sort()
    last=array[0]
    for el in array:
        if el>min_ref:
            return last
        last=el

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

def closestMax(array,max_ref):
    array=array.copy()
    array.sort()
    array.reverse()
    last=array[0]
    for el in array:
        if el<max_ref:
            return last
        last=el

def fDegrau(t_array,tamanho=1):
    y=[]
    for t in t_array:
        tmp=0
        if t>=0:
            tmp=tamanho
        y.append(tmp)
    return y

def getTimeConst(t_array,array):
    last=getLastStableEl(array)
    for i in range(len(t_array)):
        if array[i]>last*0.632:
            return t_array[i]

def closeLoop(tf,k_p=0,h_multiplier=1,D=None):
    H=control.TransferFunction([h_multiplier],[1])
    H=control.tf2io(H)
    H.name='H'
    H.set_inputs(['in'])
    H.set_outputs(['out'])

    if D==None:
        D=control.TransferFunction([k_p],[1])
    D=control.tf2io(D)
    D.name='D'
    D.set_inputs(['in'])
    D.set_outputs(['out'])

    system=control.tf2io(tf)
    system.name='system'
    system.set_inputs(['in'])
    system.set_outputs(['out'])

    #        +
    # in ---->o---->--D-->system----> out
    #         |-                 |
    #         -------H------------     
    system_MF = control.InterconnectedSystem([H,D,system] ,name='system_MF',
        connections=[
            ['H.in','system.out'],
            ['D.in','-H.out'],
            ['system.in','D.out'],
        ],
        inplist=['D.in'],
        inputs=['in'],
        outlist=['system.out','D.out','-H.out'],
        outputs=['out','control','error'])
    return system_MF


def getRootRadius(root):
    c_part=root.imag
    r_part=root.real
    radius=(r_part**2+c_part**2)**.5
    return radius

t_inicial=0
t_final=10
passo_T=0.001
T=0.1

time_array = np.arange(t_inicial, t_final,passo_T)
tfs = control.TransferFunction([1500, 16000], [750, 1500, 16000])
time_array_disc = np.arange(t_inicial, t_final,T)
tfz = control.c2d(tfs, T)
print(tfz)
print('polos:',tfz.pole())
print()

roots,gains=control.rlocus(tfz,xlim=[-20,2], kvect=np.arange(0,20,0.001), ylim=[-5,5],Plot=False)
last_gain=None
for i in range(len(roots)):
    gain=gains[i]
    root1=roots[i][0]
    root2=roots[i][1]
    radius1=getRootRadius(root1)
    radius2=getRootRadius(root2)
    if radius1<=1 and radius2<=1:
        last_gain=gain
    else:
        break

print('Hard tunning - last gain:',last_gain)
radius1=0
radius2=0
gain=last_gain
itr_count=0
bigger_than_one=False
limit_print=33
itr_max=200000
while radius1 != 1 or radius2 != 1:
    last_gain=gain
    if bigger_than_one:
        gain+=0.000000001
    else:
        gain+=((1-radius1)**.5+(1-radius2)**.5)/2
    roots,_=control.rlocus(tfz,kvect=[gain],Plot=False)
    root1=roots[0][0]
    root2=roots[0][1]
    radius1=getRootRadius(root1)
    radius2=getRootRadius(root2)
    if radius1>1 or radius2>1:
        if not bigger_than_one:
            bigger_than_one=True
            gain=last_gain
            limit_print=33333
        else:
            gain=last_gain
            root = control.rlocus(tfz,kvect=[gain],Plot=False)[0][0]
            radius1=getRootRadius(root[0])
            radius2=getRootRadius(root[1])
            break
    if itr_count+1>itr_max and bigger_than_one:
        gain=last_gain
        root = control.rlocus(tfz,kvect=[gain],Plot=False)[0][0]
        radius1=getRootRadius(root[0])
        radius2=getRootRadius(root[1])
        break
    itr_count+=1
    if itr_count%limit_print==0:
        print('Iter: {} - Radius: {} and {}'.format(itr_count,radius1,radius2))
print('Fine tunning - Kp={} gives radius={} and {} - after {} iters'.format(gain,radius1,radius2,itr_count))


print()
root,gain = control.rlocus(tfz,kvect=[9.679199999016548],Plot=False)
print('O ganho:',gain[0],'gera um polo dentro do raio, raio:',getRootRadius(root[0][0]))
system_mf=closeLoop(tfz,k_p=gain[0])
transad = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(transad[0],transad[1][0],'r',label='Resposta ao degrau unitário')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()
plt.step(transad[0],transad[1][2],'g',label='Erro')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()
plt.step(transad[0],transad[1][1],'b',label='Ação de controle')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()




roots,gains=control.rlocus(tfz,xlim=[-2.5,2.5], ylim=[-2,2])
plt.legend(loc='best')
plt.show()
print()
roots,gains=control.rlocus(tfz,Plot=False)
real_kps=[]
imag_kps=[]
for i in range(len(roots)):
    gain=gains[i]
    root1=roots[i][0]
    root2=roots[i][1]
    if root1.imag!=0 and root2.imag!=0:
        imag_kps.append(gain)
    else:
        real_kps.append(gain)

step=0.00001


if len(real_kps)<=1 or len(imag_kps)<=1:
    if len(imag_kps)<=1:
        print('Não há intervalo complexo')
    if len(real_kps)<=1:
        print('Não há intervalo real')
else:
    min_imag_kp=min(imag_kps)
    min_real_kp=min(real_kps)
    if min_real_kp < min_imag_kp:
        start=min_real_kp
        end=closestMin(imag_kps,min_real_kp)
        real_start=True
    else:
        start=min_imag_kp
        end=closestMin(real_kps,min_imag_kp)
        real_start=False

    threshold_start=start
    for kp in np.arange(start,end,step):
        root,gain = control.rlocus(tfz,kvect=[kp],Plot=False)
        gain=gain[0]
        root1=root[0][0]
        root2=root[0][1]
        is_imag=root1.imag!=0 and root2.imag!=0
        if (real_start and is_imag) or (not real_start and not is_imag):
            break
        threshold_start=kp

    if real_start:
        print('Intervalo imaginário começa em:',threshold_start)
    else:
        print('Intervalo real começa em:',threshold_start)

    max_imag_kp=max(imag_kps)
    max_real_kp=max(real_kps)
    if not(max_imag_kp < min_real_kp and not real_start) and not(max_real_kp < min_imag_kp and real_start):
        if max_real_kp < max_imag_kp:
            start=max_real_kp
            end=closestMax(imag_kps,max_real_kp)
            real_start=True
        else:
            start=max_imag_kp
            end=closestMax(real_kps,max_imag_kp)
            real_start=False

        threshold_end=start
        for kp in np.arange(start,end,step):
            root,gain = control.rlocus(tfz,kvect=[kp],Plot=False)
            gain=gain[0]
            root1=root[0][0]
            root2=root[0][1]
            is_imag=root1.imag!=0 and root2.imag!=0
            if (real_start and is_imag) or (not real_start and not is_imag):
                break
            threshold_end=kp

        if threshold_end>threshold_start:
            if real_start:
                print('Intervalo real termina em:',threshold_end)
            else:
                print('Intervalo imaginário termina em:',threshold_end)


print()
root,gain = control.rlocus(tfz,kvect=[8.45762],Plot=False)
print('O ganho:',gain[0],'gera um polo dentro do raio, raio:',getRootRadius(root[0][0]))
system_mf=closeLoop(tfz,k_p=gain[0])
transad = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(transad[0],transad[1][2],'g',label='Erro')
plt.step(transad[0],transad[1][1],'b',label='Ação de controle')
plt.step(transad[0],transad[1][0],'r',label='Resposta ao degrau unitário')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()

eu = control.forced_response(tfz,time_array_disc,fDegrau(time_array_disc))
print('Constante de tempo:',getTimeConst(eu[0],eu[1]))


print()
dahllin=control.TransferFunction([2.31461,-3.77049969,1.894971207],[1,-0.6722,0.111953,-0.632121,0.192368],T)
print(dahllin)
system_mf=closeLoop(tfz,D=dahllin)
tuca = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(tuca[0],tuca[1][2],'g',label='Erro')
plt.step(tuca[0],tuca[1][1],'b',label='Ação de controle')
plt.step(tuca[0],tuca[1][0],'r',label='Resposta ao degrau unitário')
plt.title('dahllin')
plt.legend(loc='best')
plt.show()