#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import numpy as np
import getopt
import sys
from scipy import signal,integrate
from itertools import zip_longest
import matplotlib.pyplot as plt

def getParametersMMA(){
    k=16000
    m=750
    b=1500
    return m,k,b 
}

def sqrt(x){
    return x**.5
}

def isNaN(x){
    return x != x
}

def fTransferenciaMMA_t(t){
    e=math.e
    sin=math.sin
    cos=math.cos
    y=(2*e**-t*(29*sin(sqrt(61/3)*t)+sqrt(183)*cos(sqrt(61/3)*t)))/sqrt(183)
    return y
}


def diagonalizar(M){
    size=M.shape[0]
    auto_valores,auto_vetores = np.linalg.eig(arr)
    auto_valores=auto_valores.tolist()
    auto_valores.sort()
    diag=np.identity(size)*np.array(auto_valores)
    return diag
}

def fTransferenciaMMA_s(s,diagonalizar=False){
    m,k,b=getParametersMMA()
    A=np.array([[0,1],[-k/m,-b/m]])
    B=np.array([[0,0],[k/m,b/m]])
    C=np.array([1,0])
    D=np.array([0,0])
    if diagonalizar{
        A=diagonalizar(A)
    }
    return (np.linalg.inv(C.dot((s*np.identity(2)).subtract(A))).dot(B)).add(D)
}

def fTransferenciaMMA_z(z,T,diagonalizar=False){
    m,k,b=getParametersMMA()
    A=np.array([[0,1],[((T*(b-k*T))/m)-1,2-(b*T/m)]])
    B=np.array([[0,0],[T*(k*T-b)/m,b*T/m]])
    C=np.array([1,0])
    D=np.array([0,0])
    if diagonalizar{
        A=diagonalizar(A)
    }
    return (np.linalg.inv(C.dot((z*np.identity(2)).subtract(A))).dot(B)).add(D)
}

def fImpulso(t,tamanho=1){
    y=0
    if t==0{
        y=tamanho
    }
    return y
}

def fDegrau(t,tamanho=1){
    y=0
    if t>=0{
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

def fRampa(t){
    y=t
    return y
}

def fSenoide(t){
    y=np.sin(t)
    return y
}

def removeStrPrefix(text, prefix){
    if text.startswith(prefix){
        return text[len(prefix):]
    }
    return text 
}

def atrasarSinal(a,delta){
    return np.concatenate([a[delta:],[0]*delta])
}

def somaSinais(a,b){
    y=[float(i)+float(j)  for i,j in zip_longest(a,b,fillvalue=0)]
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

def iteracaoEDOEstadoContinuo(z,t,u,u_linha=None){
    m,k,b=getParametersMMA()
    x1 = z[0]
    x2 = z[1]
    dx1dt = x2
    if u_linha is None {
        dx2dt = (-b/m)*x2-(k/m)*x1+u/m # metodo 2
    }else{
        dx2dt = (-b/m)*x2-(k/m)*x1+(b/m)*u_linha+(k/m)*u # metodo 1
    }
    dzdt = [dx1dt,dx2dt]
    return dzdt
}

def resolveEDOEstadoContinuo(t_array,array_entrada,metodo1=False){
    _,k,b=getParametersMMA()

    x1=np.empty_like(t_array)
    x2=np.empty_like(t_array)
    z0=[0,0]
    x1[0] = z0[0]
    x2[0] = z0[1]
    y=[]
    u=array_entrada

    for i in range(0,len(t_array)-1){
        tspan = [t_array[i],t_array[i+1]] # Resolve EDO de 2 em 2 passos
        extra_args= None if not metodo1 else (u[i]-u[i+1])/(t_array[i+1]-t_array[i])
        z = integrate.odeint(iteracaoEDOEstadoContinuo,z0,tspan,args=(u[i],extra_args)) # Soluciona a EDO
        # Atualiza valor de x1 e x2 a cada iteracao
        x1[i] = z[1][0]
        x2[i] = z[1][1]
        if(not metodo1){
            y.append(b*x2[i]+k*x1[i])
        }else{
            y.append(x1[i])
        }
        z0 = z[1] # proxima condicao inicial
    }
    return y,x1,x2
}

def autoVetores(A){
    auto_valores, v = np.linalg.eig(A)
    v_inv=np.linalg.inv(v)
    return v, v_inv
}

def resolveEDOEstadoDiscreto(t_array,array_entrada,diagonalizar=False){
    m,k,b=getParametersMMA()

    x1=np.empty_like(t_array)
    x2=np.empty_like(t_array)
    z0=[0,0]
    x1[0] = z0[0]
    x2[0] = z0[1]
    u=array_entrada
    T=t_array[1]-t_array[0]
    y=[]

    if(diagonalizar){
        A=np.array([[0,1],[(-1+((b*T)/m)-((k*(T**2))/m)),(2-((b*T)/m))]])
        B=np.array([0,(T**2)/m])
        C=np.array([k-(b/T),b/T])
        v, v_inv = autoVetores(A)
        A=(v_inv.dot(A)).dot(v)
        B=v_inv.dot(B)
        C=C.dot(v)
        for i in range(1,len(t_array)){
            x1[i] = A[0][0]*x1[i-1]+B[0]*u[i-1]#+A[0][1]*x2[i-1]
            x2[i] = A[1][1]*x2[i-1]+B[1]*u[i-1]#+A[1][0]*x1[i-1]
            y.append(C[0]*x1[i]+C[1]*x2[i])
        }
    }else{
        for i in range(1,len(t_array)){
            x1[i] = x2[i-1]
            x2[i] = (-1+(b*T)/m - (k*T**2)/m)*x1[i-1] + (2-(b*T)/m)*x2[i-1] + (T**2/m)*u[i-1]
            y.append((k-b/T)*x1[i-1]+(b/T)*x2[i-1])
        }
    }

    return y,x1,x2
}

def plotX1eX2(t_array,X1,X2,t_final,passo_T,nome='',discrete=False,fix_size=False,use_stem=False){
    padding_x1='    '
    unit_x2='m/s'
    if discrete {
        padding_x1=' '
        unit_x2='m'
    }
    plotArrayTemporal(t_array,X1,t_final,passo_T,label='Saída X1 |m|{}- {}'.format(padding_x1,nome),discrete=discrete,fix_size=fix_size,use_stem=use_stem)
    plotArrayTemporal(t_array,X2,t_final,passo_T,label='Saída X2 |{}| - {}'.format(unit_x2,nome),discrete=discrete,fix_size=fix_size,use_stem=use_stem)
}


def execRespostaImpulso(t_inicial,t_final,passo_T){
    t_array=[]
    func_trans_array=[]
    impulso_array=[]

    # a função fTransferenciaMMA_t não suporta valores negativos
    if t_inicial<0{ 
        t_final+=-t_inicial
        t_inicial+=-t_inicial
    }

    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        func_trans_array.append(fTransferenciaMMA_t(t))
        impulso_array.append(fImpulso(t))
    }

    res_imp = convolucao(func_trans_array, impulso_array,t_final,passo_T)

    plotArrayTemporal(t_array,res_imp,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.legend(['Resposta ao impulso'])
    plt.title('Resposta ao impulso')
    plt.show()
}

def execRespostaImpulsos(t_inicial,t_final,passo_T){
    t_array=[]
    func_trans_array=[]
    impulso_tam_2_array=[]
    impulso_tam_3_array=[]

    # a função fTransferenciaMMA_t não suporta valores negativos
    if t_inicial<0{ 
        t_final+=-t_inicial
        t_inicial+=-t_inicial
    }

    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        func_trans_array.append(fTransferenciaMMA_t(t))
        impulso_tam_2_array.append(fImpulso(t,tamanho=2))
        impulso_tam_3_array.append(fImpulso(t,tamanho=3))
    }

    res_imp_2 = convolucao(func_trans_array, impulso_tam_2_array,t_final,passo_T)
    res_imp_3 = convolucao(func_trans_array, impulso_tam_3_array,t_final,passo_T)

    plotArrayTemporal(t_array,res_imp_2,t_final,passo_T)
    plotArrayTemporal(t_array,res_imp_3,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.legend(['Resposta ao impulso tamanho 2','Resposta ao impulso tamanho 3'])
    plt.title('Resposta aos impulsos')
    plt.show()
}

def execProvaInvarianciaTemp(t_inicial,t_final,passo_T,atraso){
    t_array=[]
    impulso_array=[]
    func_trans_array=[]
    impulso_atrasado_array=[]

    # a função fTransferenciaMMA_t não suporta valores negativos
    if t_inicial<0{ 
        t_final+=-t_inicial
        t_inicial+=-t_inicial
    }

    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        func_trans_array.append(fTransferenciaMMA_t(t))
        impulso_array.append(fImpulso(t))
        impulso_atrasado_array.append(fImpulso(t-atraso))
    }

    res_imp = convolucao(func_trans_array, impulso_array,t_final,passo_T)
    res_imp_atrasado = convolucao(func_trans_array, impulso_atrasado_array,t_final,passo_T)
    res_imp_atrasado_compensado = atrasarSinal(res_imp_atrasado,int(atraso/passo_T))

    plotArrayTemporal(t_array,res_imp,t_final,passo_T)
    plotArrayTemporal(t_array,res_imp_atrasado,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.legend(['Resposta ao impulso','Resposta ao impulso atrasado'])
    plt.title('Prova Invariancia Temporal (1/2)')
    plt.show()

    plt.plot(t_array,res_imp)
    plt.plot(t_array,res_imp_atrasado_compensado)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.legend(['Resposta ao impulso','Resposta ao impulso atrasado Compensado'])
    plt.title('Prova Invariancia Temporal (2/2)')
    plt.show()
}

def execProvaLinearidade(t_inicial,t_final,passo_T){
    t_array=[]
    func_trans_array=[]
    impulso_tam_2_array=[]
    impulso_tam_3_array=[]
    impulso_tam_5_array=[]

    # a função fTransferenciaMMA_t não suporta valores negativos
    if t_inicial<0{ 
        t_final+=-t_inicial
        t_inicial+=-t_inicial
    }

    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        func_trans_array.append(fTransferenciaMMA_t(t))
        impulso_tam_2_array.append(fImpulso(t,tamanho=2))
        impulso_tam_3_array.append(fImpulso(t,tamanho=3))
        impulso_tam_5_array.append(fImpulso(t,tamanho=5))
    }

    res_imp_2 = convolucao(func_trans_array, impulso_tam_2_array,t_final,passo_T)
    res_imp_3 = convolucao(func_trans_array, impulso_tam_3_array,t_final,passo_T)
    res_imp_5 = convolucao(func_trans_array, impulso_tam_5_array,t_final,passo_T)
    res_soma_imp2_imp3=somaSinais(res_imp_2,res_imp_3)

    plotArrayTemporal(t_array,res_imp_5,t_final,passo_T)
    plotArrayTemporal(t_array,res_soma_imp2_imp3,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.legend(['Resposta ao impulso de tamanho 5','Resposta ao impulso de tamanho 2 + tamanho 3'])
    plt.title('Prova Linearidade')
    plt.show()
}

def execTendenciaEstabilidade(t_inicial,t_final,passo_T){
    t_array=[]
    func_trans_array=[]
    degrau_array=[]

    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        func_trans_array.append(fTransferenciaMMA_t(t))
        degrau_array.append(fDegrau(t))
    }

    res_degrau = convolucao(func_trans_array, degrau_array,t_final,passo_T)

    plotArrayTemporal(t_array,res_degrau,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.legend(['Resposta ao degrau'])
    plt.title('Prova Estabilidade')
    plt.show()
}

def execConvolucaoTester(t_inicial,t_final,passo_T){
    t_array=[]
    pulso_array=[]
    pulso2_array=[]

    # a convolução entre dois pulsos é melhor vizualizada em um intervalo -t t
    if t_inicial>=0{ 
        tamanho=(t_final-t_inicial)/2
        t_final-=tamanho
        t_inicial-=tamanho
    }

    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        pulso_array.append(fPulso(t,-2,2,tamanho=2))
        pulso2_array.append(fPulso(t,-1,1,tamanho=2))
    }

    conv_pulso_pulso = convolucao(pulso_array, pulso_array,t_inicial,passo_T)
    conv_pulso_pulso2 = convolucao(pulso_array, pulso2_array,t_inicial,passo_T)

    plotArrayTemporal(t_array,pulso_array,t_final,passo_T)
    plotArrayTemporal(t_array,conv_pulso_pulso,t_final,passo_T)
    plt.legend(['Pulso','Convolução'],loc='best')
    plt.title('Convolução entre dois pulsos iguais (1/2)')
    plt.show()

    plotArrayTemporal(t_array,pulso_array,t_final,passo_T)
    plotArrayTemporal(t_array,pulso2_array,t_final,passo_T)
    plotArrayTemporal(t_array,conv_pulso_pulso2,t_final,passo_T)
    plt.legend(['Pulso 1','Pulso 2','Convolução'],loc='best')
    plt.title('Convolução entre dois pulsos diferentes (2/2)')
    plt.show()
}

def execComportamentoEspacoDeEstadosContinuo(t_inicial,t_final,passo_T,metodo1=True){
    t_array=[]
    imp_array=[]
    degrau_array=[]
    rampa_array=[]
    senoide_array=[]
    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        imp_array.append(fImpulso(t))
        degrau_array.append(fDegrau(t))
        rampa_array.append(fRampa(t))
        senoide_array.append(fSenoide(t))
    }

    resp_imp,x1_imp,x2_imp=resolveEDOEstadoContinuo(t_array,imp_array,metodo1=metodo1)
    resp_degrau,x1_deg,x2_deg=resolveEDOEstadoContinuo(t_array,degrau_array,metodo1=metodo1)
    resp_rampa,x1_ramp,x2_ramp=resolveEDOEstadoContinuo(t_array,rampa_array,metodo1=metodo1)
    resp_senoide,x1_sen,x2_sen=resolveEDOEstadoContinuo(t_array,senoide_array,metodo1=metodo1)
    

    plotX1eX2(t_array,x1_imp,x2_imp,t_final,passo_T,nome='Impulso',discrete=False)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Continuo (1/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_imp,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Continuo (2/8)')
    plt.legend(['Resposta ao Impulso'],loc='best')
    plt.show()


    plotX1eX2(t_array,x1_deg,x2_deg,t_final,passo_T,nome='Degrau',discrete=False,fix_size=True)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Continuo (3/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_degrau,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Continuo (4/8)')
    plt.legend(['Resposta ao Degrau'],loc='best')
    plt.show()


    plotX1eX2(t_array,x1_ramp,x2_ramp,t_final,passo_T,nome='Rampa',discrete=False,fix_size=True)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Continuo (5/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_rampa,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Continuo (6/8)')
    plt.legend(['Resposta à Rampa'],loc='best')
    plt.show()

    plotX1eX2(t_array,x1_sen,x2_sen,t_final,passo_T,nome='Senoide',discrete=False,fix_size=True)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Continuo (7/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_senoide,t_final,passo_T)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Continuo (8/8)')
    plt.legend(['Resposta à Senoide'],loc='best')
    plt.show()
    
}

def execComportamentoEspacoDeEstadosDiscreto(t_inicial,t_final,passo_T,diagonalizar,use_stem){
    
    t_array=[]
    imp_array=[]
    degrau_array=[]
    rampa_array=[]
    senoide_array=[]
    for t in np.arange(t_inicial, t_final, passo_T){
        t_array.append(t)
        imp_array.append(fImpulso(t))
        degrau_array.append(fDegrau(t))
        rampa_array.append(fRampa(t))
        senoide_array.append(fSenoide(t))
    }

    resp_imp,x1_imp,x2_imp=resolveEDOEstadoDiscreto(t_array,imp_array,diagonalizar)
    resp_degrau,x1_deg,x2_deg=resolveEDOEstadoDiscreto(t_array,degrau_array,diagonalizar)
    resp_rampa,x1_ramp,x2_ramp=resolveEDOEstadoDiscreto(t_array,rampa_array,diagonalizar)
    resp_senoide,x1_sen,x2_sen=resolveEDOEstadoDiscreto(t_array,senoide_array,diagonalizar)
    
    plotX1eX2(t_array,x1_imp,x2_imp,t_final,passo_T,nome='Impulso',discrete=True,use_stem=use_stem)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Discreto (1/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_imp,t_final,passo_T,discrete=True,use_stem=use_stem)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Discreto (2/8)')
    plt.legend(['Resposta ao Impulso'],loc='best')
    plt.show()

    plotX1eX2(t_array,x1_deg,x2_deg,t_final,passo_T,nome='Degrau',discrete=True,use_stem=use_stem)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Discreto (3/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_degrau,t_final,passo_T,discrete=True,use_stem=use_stem)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Discreto (4/8)')
    plt.legend(['Resposta ao Degrau'],loc='best')
    plt.show()

    plotX1eX2(t_array,x1_ramp,x2_ramp,t_final,passo_T,nome='Rampa',discrete=True,use_stem=use_stem)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Discreto (5/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_rampa,t_final,passo_T,discrete=True,use_stem=use_stem)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Discreto (6/8)')
    plt.legend(['Resposta à Rampa'],loc='best')
    plt.show()

    plotX1eX2(t_array,x1_sen,x2_sen,t_final,passo_T,nome='Senoide',discrete=True,use_stem=use_stem)
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Variáveis de Espaço de Estados Discreto (7/8)')
    plt.legend(loc='best')
    plt.show()

    plotArrayTemporal(t_array,resp_senoide,t_final,passo_T,discrete=True,use_stem=use_stem)
    plt.ylabel('altura (m)')
    plt.xlabel('tempo (s)')
    plt.title('Comportamento Espaço de Estados Discreto (8/8)')
    plt.legend(['Resposta à Senoide'],loc='best')
    plt.show()

    return None
}



def main(argv){
    HELP_STR=r'Pytho{N}.py MMA-CDSD-Simulator.py [-s|--tInicial <valor>] [-e|--tFinal <valor>] [-T|--tPasso <valor>] [-d|--diagonalizar] [-a|--atraso <valor>] [-m|--metodo <valor>] [--stem] [--all] [--respostaImpulso] [--respostaImpulsos] [--provaInvarianciaTemporal] [--provaLinearidade] [--tendenciaEstabilidade] [--espacoEstadosContinuo] [--espacoEstadosDiscreto] [--testeConvolucao]'
    
    start=0
    end=8
    step=0.01
    atraso=1
    diagonalizar=False
    run_all=False
    modules_to_run=[]
    metodo1=True
    use_stem=False
    modules=["respostaImpulso","respostaImpulsos","provaInvarianciaTemporal","provaLinearidade","tendenciaEstabilidade","espacoEstadosContinuo","espacoEstadosDiscreto","testeConvolucao"]
    try{
        opts, args = getopt.getopt(argv,"hs:e:T:a:dm:",["diagonalizar","atraso=","tInicial=","tFinal","tPasso","all","metodo=", "stem"]+modules)
    }except getopt.GetoptError{
        print (HELP_STR)
        sys.exit(2)
    }
    for opt, arg in opts{
        if opt == '-h'{
            print (HELP_STR)
            sys.exit()
        }elif opt in ("-s","--tInicial"){
            start=float(arg)
        }elif opt in ("-e","--tFinal"){
            end=float(arg)
        }elif opt in ("-T","--tPasso"){
            step=float(arg)
        }elif opt in ("-d","--diagonalizar"){
            diagonalizar=True
        }elif opt in ("-a","--atraso"){
            atraso=float(arg)
        }elif opt in ("-m","--metodo"){
            if(int(arg[1:])==2){
                metodo1=False
            }
        }elif opt == "--stem"{
            use_stem=True
        }elif opt == "--all"{
            run_all=True
        }else{
            modules_to_run.append(opt)
        }
    }
    if (run_all){
        modules_to_run=modules
    }
    for module in modules_to_run{
        module=removeStrPrefix(module,'--')
        if module == "respostaImpulso"{
            execRespostaImpulso(start,end,step)
        }elif module == "respostaImpulsos"{
            execRespostaImpulsos(start,end,step)
        }elif module == "provaInvarianciaTemporal"{
            execProvaInvarianciaTemp(start,end,step,atraso)
        }elif module == "provaLinearidade"{
            execProvaLinearidade(start,end,step)
        }elif module == "tendenciaEstabilidade"{
            execTendenciaEstabilidade(start,end,step)
        }elif module == "espacoEstadosContinuo"{
            execComportamentoEspacoDeEstadosContinuo(start,end,step,metodo1)
        }elif module == "espacoEstadosDiscreto"{
            execComportamentoEspacoDeEstadosDiscreto(start,end,step,diagonalizar,use_stem=use_stem)
        }elif module == "testeConvolucao"{
            execConvolucaoTester(start,end,step)
        }else{
            print("Argumento desconhecido {}".format(module))
            print (HELP_STR)
            sys.exit(2)
        }
    }
}
    
if __name__ == "__main__"{
    main(sys.argv[1:])
}
