#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import control
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal as sgn

t_inicial=0
t_final=2
passo_T=0.01

def funcTransfer(num, den){
    return control.TransferFunction(num,den)
}

def fDegrau(tamanho=1){
    y=[]
    for t in np.arange(t_inicial, t_final, passo_T){
        tmp=0
        if t>=0{
            tmp=tamanho
        }
        y.append(tmp)
    }
    return y
}

def main(argv){
    funcao_transferencia = funcTransfer(np.array([1.]), np.array([0.05,1.]))
    time_array = np.arange(t_inicial, t_final,passo_T)


    resposta_degrau = control.forced_response(funcao_transferencia,time_array,fDegrau(5))
    plt.plot(resposta_degrau[0], fDegrau(5))
    plt.plot(resposta_degrau[0], resposta_degrau[2])

    print(resposta_degrau[2][5])

    plt.legend(['u(t)', 'y(t)'],loc='best')
    plt.xlabel('tempo (s)')
    plt.title('Transformada Laplace - Velocidade Angular x Torque aplicado | J=1 Nm/s² e b=20 NM/s')
    plt.show()

    # residuo
    # polos
    # Coefficients of the direct polynomial term. ???
    fracoes_parciais = sgn.residue(np.array([1.]), np.array([0.05,1.]), tol=0.001, rtype='avg') # 20/(s+20)

    resposta_degrau = control.forced_response(funcTransfer(np.array(fracoes_parciais[0]), [1]+(np.array(fracoes_parciais[1])*-1).tolist()),T=np.arange(t_inicial, t_final, passo_T),U=fDegrau(5))

    plt.plot(resposta_degrau[0], fDegrau(5))
    plt.plot(resposta_degrau[0], resposta_degrau[2])

    plt.legend(['u(t)', 'y(t)'],loc='best')
    plt.xlabel('tempo (s)')
    plt.title('Analítico - Velocidade Angular x Torque aplicado | J=1 Nm/s² e b=20 NM/s')
    plt.show()


    tf = ([0.01], [1.,-0.8], 0.01)
    tuca = sgn.dlsim(tf, fDegrau(5),t=time_array)

    plt.scatter(tuca[0], fDegrau(5),s=4)
    plt.scatter(tuca[0], tuca[1],s=4)
    plt.legend(['u(t)', 'y(t)'],loc='best')
    plt.xlabel('tempo (s)')
    plt.title('Transformada Z - Velocidade Angular x Torque aplicado | J=1 Nm/s² e b=20 NM/s')
    plt.show()
}

if __name__ == "__main__"{
    main(sys.argv[1:])
}