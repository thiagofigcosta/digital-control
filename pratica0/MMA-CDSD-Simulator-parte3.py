import cmath
import control
import matplotlib.pyplot as plt
import numpy as np
import math

def getFuncaoTransferenciaPulsada(s, T=0.1):
    e=math.e
    z=e**(s*T)
    return (0.2731*z-0.08311)/(z**2-1.629*z+0.8187)

def getFuncaoTransferenciaContinua(s):
    return ((1500*s+16000)/(750*(s**2)+1500*s+16000))

def evalArray(func, array):
    out=[]
    for el in array:
        out.append(func(el))
    return out

def fDegrau(array, tamanho=1):
    y=[]
    for t in array:
        tmp=0
        if t>=0:
            tmp=tamanho
        y.append(tmp)
    return y
    
t_inicial=0
t_final=5
passo_T=0.1

time_array_continuo = np.arange(t_inicial, t_final,passo_T)

time_array_discreto = np.arange(t_inicial, t_final,passo_T)


y = evalArray(getFuncaoTransferenciaContinua, time_array_continuo)
z = evalArray(getFuncaoTransferenciaPulsada, time_array_continuo)

plt.plot(time_array_continuo,y,label='H(s)')
plt.plot(time_array_continuo,z,label='G(s)')
plt.legend(loc='best')
plt.title('Compara Comportamento Frequência H(s) e G(s)')
plt.show()


tfs = control.TransferFunction([1500, 16000], [750, 1500, 16000])
tfz = control.c2d(tfs, 0.01)
tfz2 = control.c2d(tfs, 0.1)

print("Função Transferência Pulsada - T=0.01")
print(tfz)
print("")
print("----------------------------------")
print("")
print("Função Transferência Pulsada - T=0.1")
print(tfz2)

resposta_degrau = control.forced_response(tfs,time_array_continuo,fDegrau(time_array_continuo))
resposta_degrau2 = control.forced_response(tfz2,time_array_discreto,fDegrau(time_array_discreto))

plt.plot(resposta_degrau[0],resposta_degrau[1],linewidth=5,label='H(s)')
plt.step(resposta_degrau2[0],resposta_degrau2[1],label='G(z)')
plt.legend(loc='best')
plt.title('Comportamento Domínio Tempo - H(s) e G(z)')
plt.show()

k_array=[1,3,5,10]

for k_p in k_array:

    tfs_fechada = (k_p*tfs)/(1+k_p*tfs)#control.TransferFunction([1500, 16000], [750, 3000, 32000])
    resposta_degrau3 = control.forced_response(tfs_fechada,time_array_discreto,fDegrau(time_array_discreto))
    print("---------------------------------")
    print("Analógico")
    plt.plot(resposta_degrau[0],resposta_degrau[1],linewidth=5,label='Malha Aberta')
    plt.plot(resposta_degrau3[0],resposta_degrau3[1],label='Malha Fechada (Kp={})'.format(k_p))
    plt.legend(loc='best')
    plt.title('Analógico - Resposta ao Degrau')
    plt.show()


    tfz2_fechada = (k_p*tfz2)/(1+k_p*tfz2)
    resposta_degrau4 = control.forced_response(tfz2_fechada,time_array_discreto,fDegrau(time_array_discreto))
    print("---------------------------------")
    print("Digital")
    plt.step(resposta_degrau2[0],resposta_degrau2[1],label='Malha Aberta')
    plt.step(resposta_degrau4[0],resposta_degrau4[1],label='Malha Fechada (Kp={})'.format(k_p))
    plt.legend(loc='best')
    plt.title('Digital - Resposta ao Degrau')
    plt.show()


    print("---------------------------------")
    print("Analógico/Digital")
    plt.plot(resposta_degrau[0],resposta_degrau[1],linewidth=5,label='Malha Aberta - Analógico')
    plt.step(resposta_degrau2[0],resposta_degrau2[1],label='Malha Aberta - Discreto')
    plt.plot(resposta_degrau3[0],resposta_degrau3[1],label='Malha Fechada - Analógico')
    plt.step(resposta_degrau4[0],resposta_degrau4[1],label='Malha Fechada - Discreto')
    plt.legend(loc='best')
    plt.title('Resposta ao Degrau')
    plt.show()