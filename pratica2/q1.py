#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import signal as sgn
import matplotlib.pyplot as plt

def avaliarPolosEspacoEstados(R,L,C,s):
    return C*L*s*s+R*C*s+1

def matrizA(R,L,C):
    A=np.array([[0,1],[-1/(L*C),-R/L]])
    return A

def matrizB(R,L,C):
    B=np.array([[0],[1/(C*L)]])
    return B

def matrizC(R,L,C):
    C=np.array([1,0])
    return C

def matrizD(R,L,C):
    D=np.array([0])
    return D

def fDegrau(t_array,tamanho=1):
    y=[]
    for t in t_array:
        tmp=0
        if t>=0:
            tmp=tamanho
        y.append(tmp)
    return y

def autovalores(A):
    auto_valores,auto_vetores = np.linalg.eig(A)
    return auto_valores
    

def verificaPolos(R,L,C):
    A=matrizA(R,L,C)
    auto=autovalores(A)
    print('Autovalores de A: {}'.format(auto))
    print('Verficando se autovalores zeram o denominador da func de transferencia em epaço de estados...')
    for autoval in auto:
        den=avaliarPolosEspacoEstados(R,L,C,autoval)
        if np.isclose(den,0):
            print('Ok')
        else:
            print('Erro, valor diferente de zero: {}'.format(den))

# f
print('Caso 1, R=4, L=1, C=1:')
verificaPolos(R=4,L=1,C=1)
print()
print('Caso 2, R=1, L=4, C=1:')
verificaPolos(R=1,L=4,C=1)

# g
t_inicial=0
t_final=8
passo_T=0.01
time_array = np.arange(t_inicial, t_final,passo_T)
tf = (matrizA(R=4,L=1,C=1),matrizB(R=4,L=1,C=1),matrizC(R=4,L=1,C=1),matrizD(R=4,L=1,C=1))
boda = sgn.lsim(tf, U=fDegrau(time_array),T=time_array)

plt.plot(boda[0],boda[1],label='Resposta ao degrau',linewidth=4)
plt.plot(boda[0],boda[2],label='Estados')
plt.legend(loc='best')
plt.title('RLC - Espaço em estados - Caso 1, R=4, L=1, C=1')
plt.show()

fase_x=[]
fase_y=[]
for el in boda[2]:
    fase_x.append(el[0])
    fase_y.append(el[1])

plt.plot(fase_x,fase_y)
plt.title('Retrato de fase - Caso 1, R=4, L=1, C=1')
plt.show()

# h
t_inicial=0
t_final=20
passo_T=0.01
time_array = np.arange(t_inicial, t_final,passo_T)
tf = (matrizA(R=1,L=4,C=1),matrizB(R=1,L=4,C=1),matrizC(R=1,L=4,C=1),matrizD(R=1,L=4,C=1))
gereba = sgn.lsim(tf, U=fDegrau(time_array),T=time_array)

plt.plot(gereba[0],gereba[1],label='Resposta ao degrau',linewidth=4)
plt.plot(gereba[0],gereba[2],label='Estados')
plt.legend(loc='best')
plt.title('RLC - Espaço em estados - Caso 2, R=1, L=4, C=1')
plt.show()

fase_x=[]
fase_y=[]
for el in gereba[2]:
    fase_x.append(el[0])
    fase_y.append(el[1])

plt.plot(fase_x,fase_y)
plt.title('Retrato de fase - Caso 2, R=1, L=4, C=1')
plt.show()
