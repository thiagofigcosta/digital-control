#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import getopt
import sys
from scipy import signal,integrate
from itertools import zip_longest
import matplotlib.pyplot as plt

def fImpulso(t,tamanho=1){
    y=0
    if t==0{
        y=tamanho
    }
    return y
}

def fPulso(t,t1,t2,tamanho=1){
    y=0
    if t>=t1 and t<=t2{
        y=tamanho
    }
    return y
}

def convolucao(a,b,t_inicio=0,passo_T=0.1){
    start_idx=0
    if t_inicio < 0{
        start_idx=int((-t_inicio)/passo_T)
    }
    return (signal.convolve(a, b,method='direct')*passo_T).tolist()[start_idx:]
}

COLOR_ID=0
def plotArrayTemporal(t_array,y_array,t_final,passo_T,label=None,discrete=False,fix_size=False,use_stem=False){
    global COLOR_ID
    if fix_size {
        y_array=y_array[:-4]
    }
    if len(t_array)<len(y_array){
        t_atual=t_final
        while len(t_array)<len(y_array){
            t_array.append(t_atual)
        }
        t_atual+=passo_T
    }
    if len(t_array)>len(y_array){
        while len(t_array)>len(y_array){
            t_array.pop()
        }
    }
    if discrete{
        if use_stem{
            color=['b','g','r','m','y']
            current_color=color[COLOR_ID%len(color)]
            plt.stem(t_array,y_array,current_color,markerfmt='{}.'.format(current_color),label=label)
            COLOR_ID+=1
        }else{
            plt.scatter(t_array,y_array,s=.7,label=label)
        }
    }else{
        plt.plot(t_array,y_array,label=label)
    }
}

def main(argv){

    t_inicial=0
    t_final=31
    passo_T=1
    d=0.1

    t_array=[]

    imp_array=[]
    resp_imp_array=[]

    pulso_array=[]
    resp_pulso_array=[]
    y=[0]
    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)

        imp_array.append(fImpulso(t))
        resp_imp_array.append((1-d)**t)

        pulso_array.append(fPulso(t,1,8,tamanho=50))
    }

    for t in np.arange(t_inicial, t_final-1, passo_T){
        y.append((1-d)*y[t]+pulso_array[t+1])
    }

    

    resp_pulso_array = convolucao(pulso_array,resp_imp_array,t_final,passo_T)

    plt.stem(t_array,resp_imp_array,'g','g.',label='h(k)=(1-d)^k')
    plt.xlabel('tempo (dia)')
    plt.title('Resposta ao Impulso Conta Corrente para d=0.1')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array, resp_pulso_array, t_final, passo_T, discrete=True, use_stem=True)
    plt.stem(t_array,resp_pulso_array,'g','g.',label='f(k)=h(k)*u(k)')
    plt.xlabel('tempo (dia)')
    plt.title('Resposta ao Pulso de 50 reais por 7 dias | d=0.1')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array, y, t_final, passo_T, discrete=True, use_stem=True)
    plt.xlabel('tempo (dia)')
    plt.title('Resposta ao Pulso de 50 reais por 7 dias | Simulação Numérica | d=0.1')
    plt.legend(loc='best')
    plt.show()

}

if __name__ == "__main__"{
    main(sys.argv[1:])
}