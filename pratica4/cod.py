#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import signal as sgn
import matplotlib.pyplot as plt
import control
import sympy as sym

def fDegrau(t_array,tamanho=1):
    y=[]
    for t in t_array:
        tmp=0
        if t>=0:
            tmp=tamanho
        y.append(tmp)
    return y

def getRootRadius(root):
    c_part=root.imag
    r_part=root.real
    radius=(r_part**2+c_part**2)**.5
    return radius

def closestMin(array,min_ref):
    array=array.copy()
    array.sort()
    last=array[0]
    for el in array:
        if el>min_ref:
            return last
        last=el

def closestMax(array,max_ref):
    array=array.copy()
    array.sort()
    array.reverse()
    last=array[0]
    for el in array:
        if el<max_ref:
            return last
        last=el

def gz(z):
    return (0.004683*z + 0.004381)/(z**2 - 1.819*z + 0.8187)


def closeLoop(tf,k_p,h_multiplier=1):
    H=control.TransferFunction([h_multiplier],[1])
    H=control.tf2io(H)
    H.name='H'
    H.set_inputs(['in'])
    H.set_outputs(['out'])

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


# a
t_inicial=0
t_final=10
passo_T=0.001
T=0.1

time_array = np.arange(t_inicial, t_final,passo_T)
tfs = control.TransferFunction([.5],[.5,1,0])
time_array_disc = np.arange(t_inicial, t_final,T)
tfz = control.c2d(tfs, T)
print(tfz)

# b
roots,gains=control.rlocus(tfz,xlim=[-20,2], ylim=[-5,5])
plt.legend(loc='best')
plt.show()

# c
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

# d - tentativa 1
z = sym.symbols('z')
derivative_of_gz=sym.diff(gz(z),z)
print('Derivada')
print(derivative_of_gz)
solutions=sym.solve(derivative_of_gz,z)
print()
print('Pontos de partida:',solutions)
print()
print('Intervalo do ganho:')
for sol in solutions:
    print(1/gz(complex(sol)))

# d - tentativa 2
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


# e
print()
initial_search_space=np.arange(-500,50,1)
roots,gains=control.rlocus(tfz,kvect=initial_search_space,Plot=False)
not_oscl_kps=[]
oscl_kps=[]
for i in range(len(roots)):
    gain=gains[i]
    root1=roots[i][0]
    root2=roots[i][1]
    if (root1.imag!=0 and root2.imag!=0) or (root1.real<0 and root2.real<0):
        oscl_kps.append(gain)
    else:
        not_oscl_kps.append(gain)

step=0.00001


if len(not_oscl_kps)<=1 or len(oscl_kps)<=1:
    if len(oscl_kps)<=1:
        print('Não há intervalo oscilatório')
    if len(not_oscl_kps)<=1:
        print('Não há intervalo não oscilatório')
else:

    min_oscl_kp=min(oscl_kps)
    min_not_oscl_kp=min(not_oscl_kps)
    if min_not_oscl_kp < min_oscl_kp:
        start=min_not_oscl_kp
        end=closestMin(oscl_kps,min_not_oscl_kp)
        not_oscl_start=True
    else:
        start=min_oscl_kp
        end=closestMin(not_oscl_kps,min_oscl_kp)
        not_oscl_start=False

    threshold_start=start
    for kp in np.arange(end,start,-step): # revesed to go faster
        root,gain = control.rlocus(tfz,kvect=[kp],Plot=False)
        gain=gain[0]
        root1=root[0][0]
        root2=root[0][1]
        is_oscl=(root1.imag!=0 and root2.imag!=0) or (root1.real<0 and root2.real<0)
        # revesed to go faster
        if (not not_oscl_start and is_oscl) or (not_oscl_start and not is_oscl):
            break
        threshold_start=kp

    if not_oscl_start:
        print('Intervalo oscilatório começa em:',threshold_start)
    else:
        print('Intervalo não oscilatório começa em:',threshold_start)

    max_oscl_kp=max(oscl_kps)
    max_not_oscl_kp=max(not_oscl_kps)

    if max_not_oscl_kp < max_oscl_kp:
        start=max_not_oscl_kp
        end=closestMax(oscl_kps,max_not_oscl_kp)
        not_oscl_start=True
    else:
        start=max_oscl_kp
        end=closestMax(not_oscl_kps,max_oscl_kp)
        not_oscl_start=False

    threshold_end=start
    for kp in np.arange(start,end,step):
        root,gain = control.rlocus(tfz,kvect=[kp],Plot=False)
        gain=gain[0]
        root1=root[0][0]
        root2=root[0][1]
        is_oscl=root.imag!=0 or root.real<0
        if (not_oscl_start and is_oscl) or (not not_oscl_start and not is_oscl):
            break
        threshold_end=kp
    if threshold_end>threshold_start:
        if not_oscl_start:
            print('Intervalo não oscilatório termina em:',threshold_end)
        else:
            print('Intervalo oscilatório termina em:',threshold_end)


#f
print()
root,gain = control.rlocus(tfz,kvect=[0.94],Plot=False)
print('O ganho:',gain[0],'gera os polos reais diferentes:',root[0])
system_mf=closeLoop(tfz,k_p=gain[0])
joaoluzia = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(joaoluzia[0],joaoluzia[1][0],label='Resposta ao degrau unitário')
plt.step(joaoluzia[0],joaoluzia[1][1],label='Ação de controle')
plt.step(joaoluzia[0],joaoluzia[1][2],label='Erro')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()


# g
print()
root,gain = control.rlocus(tfz,kvect=[0.951447061],Plot=False)
print('O ganho:',gain[0],'gera dois polos reais similares:',root[0])
print('Devido a erros de arredondamento não foi possível encontrar valores exatamente iguais')
system_mf=closeLoop(tfz,k_p=gain[0])
metade = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(metade[0],metade[1][0],label='Resposta ao degrau unitário')
plt.step(metade[0],metade[1][1],label='Ação de controle')
plt.step(metade[0],metade[1][2],label='Erro')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()


#h
print()
root,gain = control.rlocus(tfz,kvect=[40],Plot=False)
print('O ganho:',gain[0],'gera os polos complexos:',root[0])
system_mf=closeLoop(tfz,k_p=gain[0])
monte = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(monte[0],monte[1][0],label='Resposta ao degrau unitário')
plt.step(monte[0],monte[1][1],label='Ação de controle')
plt.step(monte[0],monte[1][2],label='Erro')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()

#i
print()
root,gain = control.rlocus(tfz,kvect=[1600],Plot=False)
print('O ganho:',gain[0],'gera os polos reais e negativos:',root[0])
system_mf=closeLoop(tfz,k_p=gain[0])
pelad = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(pelad[0],pelad[1][0],label='Resposta ao degrau unitário')
plt.step(pelad[0],pelad[1][2],label='Erro')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()
plt.step(pelad[0],pelad[1][0],label='Resposta ao degrau unitário')
plt.step(pelad[0],pelad[1][2],label='Erro')
plt.step(pelad[0],pelad[1][1],label='Ação de controle')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()

#j
print()
root,gain = control.rlocus(tfz,kvect=[48],Plot=False)
print('O ganho:',gain[0],'gera um polo fora do raio, raio:',getRootRadius(root[0][0]))
system_mf=closeLoop(tfz,k_p=gain[0])
furunc = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(furunc[0],furunc[1][0],label='Resposta ao degrau unitário')
plt.step(furunc[0],furunc[1][1],label='Ação de controle')
plt.step(furunc[0],furunc[1][2],label='Erro')
plt.title('Kp: {}'.format(gain[0]))
plt.legend(loc='best')
plt.show()

#k
print()
kps0=[]
errors0=[]
degrau=fDegrau(time_array_disc)
for k_p in np.arange(0,41,1):
    system_mf=closeLoop(tfz,k_p=k_p)
    cagad = control.input_output_response(system_mf, time_array_disc, degrau)
    error_sum=0
    for i in range(len(cagad[0])):
        error_sum+=abs(degrau[i]-cagad[1][0][i])
    errors0.append(error_sum)
    kps0.append(k_p)
kps1=[]
errors1=[]
degrau=fDegrau(time_array_disc)
for k_p in np.arange(1,6,0.01):
    system_mf=closeLoop(tfz,k_p=k_p)
    cagad = control.input_output_response(system_mf, time_array_disc, degrau)
    error_sum=0
    for i in range(len(cagad[0])):
        error_sum+=abs(degrau[i]-cagad[1][0][i])
    errors1.append(error_sum)
    kps1.append(k_p)
kps2=[]
errors2=[]
degrau=fDegrau(time_array_disc)
for k_p in np.arange(4.9,6.1,0.001):
    system_mf=closeLoop(tfz,k_p=k_p)
    cagad = control.input_output_response(system_mf, time_array_disc, degrau)
    error_sum=0
    for i in range(len(cagad[0])):
        error_sum+=abs(degrau[i]-cagad[1][0][i])
    errors2.append(error_sum)
    kps2.append(k_p)
kps3=[]
errors3=[]
degrau=fDegrau(time_array_disc)
for k_p in np.arange(6,7,0.00001):
    system_mf=closeLoop(tfz,k_p=k_p)
    cagad = control.input_output_response(system_mf, time_array_disc, degrau)
    error_sum=0
    for i in range(len(cagad[0])):
        error_sum+=abs(degrau[i]-cagad[1][0][i])
    errors3.append(error_sum)
    kps3.append(k_p)


best=min(errors3)
best_kp=kps3[errors3.index(best)]

print('Melhor Kp:',best_kp)

plt.plot(kps0,errors0,label='Erro')
plt.title('Evolução do erro com o aumento de Kp')
plt.legend(loc='best')
plt.show()

plt.plot(kps1,errors1,label='Erro')
plt.title('Evolução do erro com o aumento de Kp (ampliado)')
plt.legend(loc='best')
plt.show()

plt.plot(kps2,errors2,label='Erro')
plt.title('Evolução do erro com o aumento de Kp (ampliado*2)')
plt.legend(loc='best')
plt.show()

plt.plot(kps3,errors3,label='Erro')
plt.title('Evolução do erro com o aumento de Kp (ampliado*3)')
plt.legend(loc='best')
plt.show()

system_mf=closeLoop(tfz,k_p=best_kp)
kaidbok = control.input_output_response(system_mf, time_array_disc, fDegrau(time_array_disc))
plt.step(kaidbok[0],kaidbok[1][0],label='Resposta ao degrau unitário')
plt.step(kaidbok[0],kaidbok[1][1],label='Ação de controle')
plt.step(kaidbok[0],kaidbok[1][2],label='Erro')
plt.title('Melhor Kp, Kp: {}'.format(best_kp))
plt.legend(loc='best')
plt.show()