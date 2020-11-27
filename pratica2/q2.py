#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import control
from scipy import signal as sgn
import matplotlib.pyplot as plt

taxa_natalidade=0.01
taxa_mortalidade_larva=0.1
taxa_mortalidade_casulo=0.1
taxa_mortalidade_besouro=0.9

def avaliarPolosEspacoEstados(u_l,u_c,u_b,v,s):
    return ((s**3)+u_b*(s**2)-(s**2)+u_c*v-u_c*u_l*v+u_l*v-v)

def matrizA(u_l,u_c,u_b,v):
    A=np.array([[0,0, v],[1-u_l,0,0],[0,1-u_c,1-u_b]])
    return A

def matrizB(u_l,u_c,u_b,v):
    B=np.array([[0],[0],[1]])
    return B

def matrizC(u_l,u_c,u_b,v):
    C=np.array([0,0,1])
    return C

def matrizD(u_l,u_c,u_b,v):
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
    return auto_valores,auto_vetores

def extraiVariaveis(variaveis):
    x=[]
    for i in range(len(variaveis[0])):
        x.append([])
    for line in variaveis:
        for j, column in enumerate(line):
            x[j].append(column)
    return x

    

def verificaPolos(A, u_l,u_c,u_b,v):
    auto,_=autovalores(A)
    print('Autovalores de A: {}'.format(auto))
    print('Verficando se autovalores zeram o denominador da func de transferencia em epaço de estados...')
    for autoval in auto:
        den=avaliarPolosEspacoEstados(u_l,u_c,u_b,v,autoval)
        if np.isclose(den,0):
            print('Ok')
        else:
            print('Erro, valor diferente de zero: {}'.format(den))


print('Caso 1 - taxa_mortalidade_larva={}, taxa_mortalidade_casulo={}, taxa_mortalidade_besouro={}, taxa_natalidade={}'.format(taxa_mortalidade_larva, taxa_mortalidade_casulo, taxa_mortalidade_besouro,taxa_natalidade))
A=matrizA(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade)
B=matrizB(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade)
C=matrizC(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade)
D=matrizD(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade)

verificaPolos(A, taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade)


t_inicial=0
t_final=8
passo_T=0.1
time_array = np.arange(t_inicial, t_final,passo_T)
tf = (A, B, C, D)
boda = sgn.lsim(tf, U=fDegrau(time_array),T=time_array)

x1,x2,x3 = extraiVariaveis(boda[2])


plt.scatter(boda[0],boda[1],label='Resposta ao degrau',s=14)
plt.scatter(boda[0],x3,label='X_3=B(k)=Nº Besouro',s=4)
plt.scatter(boda[0],x2,label='X_2=C(k)=Qtd Casulo',s=4)
plt.scatter(boda[0],x1,label='X_1=L(k)=Qtd Larva',s=4)
plt.legend(loc='best')
plt.title('Besouro - Espaço em estados - u_l={}, u_c={}, u_b={}, v={}'.format(taxa_mortalidade_larva, taxa_mortalidade_casulo, taxa_mortalidade_besouro,taxa_natalidade))
plt.show()

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(x1, x2, x3,s=4)
ax.set_xlabel('X1')
ax.set_ylabel('X2')
ax.set_zlabel('X3')
plt.title('Retrato de fase - u_l={}, u_c={}, u_b={}, v={}'.format(taxa_mortalidade_larva, taxa_mortalidade_casulo, taxa_mortalidade_besouro,taxa_natalidade))
plt.show()


stateSpace = control.StateSpace(
    matrizA(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade),
    matrizB(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade),
    matrizC(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade),
    matrizD(taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade)
)

# ’reachable’ - reachable canonical form
# ’observable’ - observable canonical form
# ’modal’ - modal canonical form
# canonical = control.canonical_form(stateSpace, form='modal')
# print('A_w: \n{}'.format(canonical[0].A))
# print('------')
# print('B_w: \n{}'.format(canonical[0].B))
# print('------')
# print('C_w: \n{}'.format(canonical[0].C))
# print('------')
# print('D_w: \n{}'.format(canonical[0].D))
# print('------')
print('Caso Diagonalização - taxa_mortalidade_larva={}, taxa_mortalidade_casulo={}, taxa_mortalidade_besouro={}, taxa_natalidade={}'.format(taxa_mortalidade_larva, taxa_mortalidade_casulo, taxa_mortalidade_besouro,taxa_natalidade))
print()
_, V = autovalores(A)
P=np.linalg.inv(V)

A_w=P.dot(A).dot(V)
B_w=P.dot(B)
C_w=C.dot(V)
D_w=D

print('V= \n{}'.format(V))
print('------')
print('P=V⁻¹= \n{}'.format(P))
print('------')
print('A_w=V⁻¹AV= \n{}'.format(A_w))
print('------')
print('B_w=PB= \n{}'.format(B_w))
print('------')
print('C_w=CP⁻¹= \n{}'.format(C_w))
print('------')
print('D_w=D= \n{}'.format(D_w))
print('------')

verificaPolos(P.dot(A).dot(V), taxa_mortalidade_larva,taxa_mortalidade_casulo,taxa_mortalidade_besouro,taxa_natalidade)

tf = (
        A_w,
        B_w,
        C_w,
        D_w
    )
boda = sgn.lsim(tf, U=fDegrau(time_array),T=time_array)

x1,x2,x3 = extraiVariaveis(boda[2])

plt.scatter(boda[0],boda[1],label='Resposta ao degrau',s=14)
plt.scatter(boda[0],x3,label='X_3',s=4)
plt.scatter(boda[0],x2,label='X_2',s=14)
plt.scatter(boda[0],x1,label='X_1',s=4)
plt.legend(loc='best')
plt.title('Besouro - Espaço em Estados - Canônico - u_l={}, u_c={}, u_b={}, v={}'.format(taxa_mortalidade_larva, taxa_mortalidade_casulo, taxa_mortalidade_besouro,taxa_natalidade))
plt.show()