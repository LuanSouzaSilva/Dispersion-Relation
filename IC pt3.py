import matplotlib.pyplot as plt
import numpy as np
import math
import time
import collections

#Declaracao de constantes, parametros e vetores
e = 2.718281828459045235360287
pi = np.pi
"""
w0 = int(input('Digite a frequencia natural do sistema: '))
nu = int(input('Digite o numero de passos: '))
N = int(input('Digite o numero de massas: '))
T = float(input('Digite quanto tempo durara a simulacao: '))
t = float(input('Digite o instante de tempo: '))
a = float(input('Digite a separacao entre as massas: '))
"""
w0 = 20
nu = 10000 
N = 3
T = 20
t = 18
a = 0.1
a0 = 0.1

N1 = 100

dt = T/nu
N2 = 200
udot = np.zeros((N, nu))  
u = np.zeros((N, nu))
A = np.zeros(N2)

k = np.zeros(N2)

k1 = 0
k2 = 0
k3 = 0
k4 = 0

m1 = np.zeros(N)
m2 = np.zeros(N)
m3 = np.zeros(N)
m4 = np.zeros(N)

im = 0 + 1j

t0 = time.time()
tempo = np.linspace(0, T, nu)
w = np.linspace(0, 2*w0, N1)


#Funcao que ira resolver o sistema de N equacoes usando o metodo RK4
def Solucao(n, num, W, W0, h):
    for j in range(num-1):
        for i in range(n):
            if i == 0:
                m1[i] = udot[i][j]*h
                k1 = h*(a0*np.cos(W*j*h)+(W0**2)*u[-1][j] - 2*(W0**2)*u[i][j] + (W0**2)*u[i+1][j])
                
                m2[i] = (udot[i][j] + k1/2)*h
                k2 = h*(a0*np.cos(h*W*(j + 0.5))+(W0**2)*(u[-1][j] + m1[-1]/2) - 2*(W0**2)*(u[i][j] + m1[i]/2) + (W0**2)*(u[i+1][j] + m1[i+1]/2))
                
                m3[i] = (udot[i][j] + k2/2)*h
                k3 = h*(a0*np.cos(h*W*(j + 0.5))+(W0**2)*(u[-1][j] + m2[-1]/2) - 2*(W0**2)*(u[i][j] + m2[i]/2) + (W0**2)*(u[i+1][j] + m2[i+1]/2))
                
                m4[i] = (udot[i][j] + k3)*h
                k4 = h*(a0*np.cos(h*W*(j + 1))+(W0**2)*(u[-1][j] + m3[-1]) - 2*(W0**2)*(u[i][j] + m3[i]) + (W0**2)*(u[i+1][j] + m3[i+1]))
                
                udot[i][j+1] = udot[i][j] + (k1 + k4 + 2*(k2 + k3))/6
                            
                u[i][j+1] = u[i][j] + (m1[i] + m4[i] + 2*(m2[i] + m3[i]))/6
                
            elif i + 1 == n:
                
                m1[i] = udot[i][j]*h
                k1 = h*((W0**2)*u[i-1][j] - 2*(W0**2)*u[i][j] + (W0**2)*u[0][j])
                
                m2[i] = (udot[i][j] + k1/2)*h
                k2 = h*((W0**2)*(u[i-1][j] + m1[i-1]/2) - 2*(W0**2)*(u[i][j] + m1[i]/2) + (W0**2)*(u[0][j] + m1[0]/2))
                
                m3[i] = (udot[i][j] + k2/3)*h
                k3 = h*((W0**2)*(u[i-1][j] + m2[i-1]/2) - 2*(W0**2)*(u[i][j] + m2[i]/2) + (W0**2)*(u[0][j] + m2[0]/2))
                
                m4[i] = (udot[i][j] + k3)*h
                k4 = h*((W0**2)*(u[i-1][j] + m3[i-1]) - 2*(W0**2)*(u[i][j] + m3[i]) + (W0**2)*(u[0][j] + m3[0]))
                
                udot[i][j+1] = udot[i][j] + (k1 + k4 + 2*(k2 + k3))/6
                
                u[i][j+1] = u[i][j] + (m1[i] + m4[i] + 2*(m2[i] + m3[i]))/6
                
                
            else:
                m1[i] = udot[i][j]*h
                k1 = h*((W0**2)*u[i-1][j] - 2*(W0**2)*u[i][j] + (W0**2)*u[i+1][j])
                
                m2[i] = (udot[i][j] + k1/2)*h
                k2 = h*((W0**2)*(u[i-1][j] + m1[i-1]/2) - 2*(W0**2)*(u[i][j] + m1[i]/2) + (W0**2)*(u[i+1][j] + m1[i+1]/2))
                
                m3[i] = (udot[i][j] + k2/3)*h
                k3 = h*((W0**2)*(u[i-1][j] + m2[i-1]/2) - 2*(W0**2)*(u[i][j] + m2[i]/2) + (W0**2)*(u[i+1][j] + m2[i+1]/2))
                
                m4[i] = (udot[i][j] + k3)*h
                k4 = h*((W0**2)*(u[i-1][j] + m3[i-1]) - 2*(W0**2)*(u[i][j] + m3[i]) + (W0**2)*(u[i+1][j] + m3[i+1]))
                
                udot[i][j+1] = udot[i][j] + (k1 + k4 + 2*(k2 + k3))/6
                
                u[i][j+1] = u[i][j] + (m1[i] + m4[i] + 2*(m2[i] + m3[i]))/6
                       
    return u


t = int(t/dt)

k1 = []
w1 = []
TFX = []
N2 = 200

for m in range(N1):
    A = np.zeros(N2)
    
    u = Solucao(N, nu, w[m], w0, dt)
       
    for i in range(N2):
        k[i] = pi*i/((N2-1)*a)
        for j in range(N-1):
            A[i] += u[j][t]*(e**(im*k[i]*j*a))
    for i in range(N2):
        A[i] = np.abs(A[i])
    B = zip(A, k)
    B = list(B)
    B = max(B)
    TFX.append(B[0])
    k1.append(B[1])
    w1.append(w[m])

w2 = []
k2 = []

count = collections.Counter(k1)
count = list(count)
TFX1 = []
w3 = []

for i in range(len(count)):
    w2 += [[]]
    TFX1 += [[]]
    for j in range(N1):       
        if k1[j] == count[i]:
            w2[i].append(w1[j])
            TFX1[i].append(TFX[j])

for i in range(len(count)):
    D = list(zip(TFX1[i], w2[i]))
    D = max(D)
    w3.append(D[1])
    #k2.append(count[i])

plt.figure(1)
plt.title('Gráfico da Relação de Dispersão', fontsize = 22)
plt.xlabel('Número de onda $k$', fontsize = 20)
plt.ylabel('Frequência angular $\omega$', fontsize = 20)

plt.plot(count, w3, 'o')        

x = np.linspace(0, pi/a, 1000)
y = 40*np.sin(x*a/2)
plt.plot(x, y)

#Pedaco opcional do codigo, que plota u(t)
#u = Solucao(N, nu, w[30], w0, dt)
#plt.figure(2)
#plt.title('Descolamento em função do tempo', fontsize = 22)
#plt.xlabel('Tempo', fontsize = 20)
#plt.ylabel('Deslocamento', fontsize = 20)
#for i in range(N):
#    plt.plot(tempo, u[i])


#plt.xlim(0, T)
#print(u[2])
#plt.show()


tf = time.time()
print('Tempo de processamento (minutos): ', (tf-t0)/60)
