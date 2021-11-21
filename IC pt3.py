import matplotlib.pyplot as plt
import numpy as np
import math
import time

#Declaracao de constantes, parametros e vetores
e = 2.718281828
pi = 3.14159265
"""
w0 = int(input('Digite a frequencia natural do sistema: '))
nu = int(input('Digite o numero de passos: '))
N = int(input('Digite o numero de massas: '))
T = float(input('Digite quanto tempo durara a simulacao: '))
t = float(input('Digite o instante de tempo: '))
a = float(input('Digite a separacao entre as massas: '))
"""
w0 = 20
nu = 9000#0
N = 2
T = 10
t = 2.99
a = 1 
a0 = 0.1

N1 = 100

N += 1

vgm = 0

dt = T/nu

udot = np.zeros((N, nu))  
u = np.zeros((N, nu))
tempo = np.zeros(nu)
A = np.zeros(N)

w = np.zeros(N1)
k = np.zeros(N)
vg = np.zeros(N)

k1 = 0
k2 = 0
k3 = 0
k4 = 0


im = 0 + 1j

t0 = time.time()

#Funcao que ira resolver o sistema de N equacoes usando o metodo RK4
def Solucao(n, num, W, W0, h):
    for j in range(num-1):
        for i in range(n-1):
            if i == 0:
                if j == 0:
                    k1 = a0*np.cos(W*j*h)+ (W0**2)*u[n-2][j] - 2*(W0**2)*u[i][j] + (W0**2)*u[i+1][j]
                    k2 = a0*np.cos(h*W*(j + 0.5))+ (W0**2)*(u[n-2][j] + k1*h/2) - 2*(W0**2)*(u[i][j] + k1*h/2) + (W0**2)*(u[i+1][j] + k1*h/2)
                    k3 = a0*np.cos(h*W*(j + 0.5))+ (W0**2)*(u[n-2][j] + k2*h/2) - 2*(W0**2)*(u[i][j] + k2*h/2) + (W0**2)*(u[i+1][j] + k2*h/2)
                    k4 = a0*np.cos(h*W*(j + 1))+ (W0**2)*(u[n-2][j] + k3*h) - 2*(W0**2)*(u[i][j] + k3*h) + (W0**2)*(u[i+1][j] + k3*h)
            
                
                    udot[i][j+1] = udot[i][j] + (k1 + k4 + 2*(k2 + k3))*h/6
                
                
                    u[i][j+1] = u[i][j] + udot[i][j]*h#(k11[i] + k41[i] + 2*(k21[i] + k31[i]))*h/6
                else:
                    k1 = a0*np.cos(W*j*h)+(W0**2)*u[n-2][j] - 2*(W0**2)*u[i][j] + (W0**2)*u[i+1][j]
                    k2 = a0*np.cos(h*W*(j + 0.5))+(W0**2)*(u[n-2][j] + k1*h/2) - 2*(W0**2)*(u[i][j] + k1*h/2) + (W0**2)*(u[i+1][j] + k1*h/2)
                    k3 = a0*np.cos(h*W*(j + 0.5))+(W0**2)*(u[n-2][j] + k2*h/2) - 2*(W0**2)*(u[i][j] + k2*h/2) + (W0**2)*(u[i+1][j] + k2*h/2)
                    k4 = a0*np.cos(h*W*(j + 1))+(W0**2)*(u[n-2][j] + k3*h) - 2*(W0**2)*(u[i][j] + k3*h) + (W0**2)*(u[i+1][j] + k3*h)
            
                
                    udot[i][j+1] = udot[i][j] + (k1 + k4 + 2*(k2 + k3))*h/6
                
                
                    u[i][j+1] = u[i][j] + udot[i][j]*h
            elif i + 2 == n:
                k1 = (W0**2)*u[i-1][j] - 2*(W0**2)*u[i][j] + (W0**2)*u[0][j]
                k2 = (W0**2)*(u[i-1][j] + k1*h/2) - 2*(W0**2)*(u[i][j] + k1*h/2) + (W0**2)*(u[0][j] + k1*h/2)
                k3 = (W0**2)*(u[i-1][j] + k2*h/2) - 2*(W0**2)*(u[i][j] + k2*h/2) + (W0**2)*(u[0][j] + k2*h/2)
                k4 = (W0**2)*(u[i-1][j] + k3*h) - 2*(W0**2)*(u[i][j] + k3*h) + (W0**2)*(u[0][j] + k3*h)
                
                udot[i][j+1] = udot[i][j] + (k1 + k4 + 2*(k2 + k3))*h/6
                
                u[i][j+1] = u[i][j] + udot[i][j]*h
                
                
            else:
                
                k1 = (W0**2)*u[i-1][j] - 2*(W0**2)*u[i][j] + (W0**2)*u[i+1][j]
                k2 = (W0**2)*(u[i-1][j] + k1*h/2) - 2*(W0**2)*(u[i][j] + k1*h/2) + (W0**2)*(u[i+1][j] + k1*h/2)
                k3 = (W0**2)*(u[i-1][j] + k2*h/2) - 2*(W0**2)*(u[i][j] + k2*h/2) + (W0**2)*(u[i+1][j] + k2*h/2)
                k4 = (W0**2)*(u[i-1][j] + k3*h) - 2*(W0**2)*(u[i][j] + k3*h) + (W0**2)*(u[i+1][j] + k3*h)
                
                udot[i][j+1] = udot[i][j] + (k1 + k4 + 2*(k2 + k3))*h/6
                
                u[i][j+1] = u[i][j] + udot[i][j]*h
                       
    return u

for i in range(N1-1):
    w[i+1] = w[i] + 2*w0/N1

"""
u = Solucao(N, nu, w[50], w0, dt)


t = int(t/dt)

b = 0
N2 = 100
"""
k1 = []
w1 = []
TFX = []

X = np.zeros((N-1, nu))
"""
for i in range(N-1):
    for j in range(nu):
        X[i][j] = i*a + u[i][j]
       
A = np.fft.fft(X[1])
A = 2*np.abs(A/nu)
k = np.fft.fftfreq(nu)
k *= nu/(20)
"""
A2 = np.zeros((N1, nu))

for m in range(N1):
    l = int(nu/2)
    b = 0
    u = Solucao(N, nu, w[m], w0, dt)
    for i in range(N-1):
        for j in range(nu):
            X[i][j] = i*a + u[i][j]

    A = np.fft.fft(X[1])
    A = 2*np.abs(A/nu)
    k = np.fft.fftfreq(nu)
    k *= nu/20
    
    A2[m] = A
    
k2 = []
C = []
TFx1 = []
for i in range(N1):
    c = 0
    for j in range(nu-1):
        if A2[i][j] > A2[i][j-1] and A2[i][j] > A2[i][j+1] and k[j]>= -3.15 and k[j] <= 3.15 and k[j] != 0:
            w1.append(w[i])
            k1.append(k[j])
            TFx1.append(A2[i][j])
            c += 1
    if c > 0:
        k2 += [[]]
        C.append(c)

#print(w1)
w2 = []
#k2 = []
TFx2 = [] 
    
for i in range(len(C)):
    for j in range(C[i]):
        k2[i].append(k1[j])
print(k2)
for i in range(N1):
    for j in range(len(w1)):
        if w[j] == w1[i]:
            w2.append(w[i])
            break

"""          
    for j in range(len(k1)):
        a = comprimento[i]
        while j > 0 and comprimento[j-1] > a:
            comprimento[j] = comprimento[j-1]
            i -= 1
        comprimento[j] = a
"""
#plt.plot(k1, w1, 'o')


"""
plt.xlim(0, pi/a)
plt.ylim(0, 0.003)
plt.grid()
plt.plot(k, A)
"""

"""
#Pedaco opcional do codigo, que plota u(t)

#contagem do tempo para o plot
for j in range(nu):
    tempo[j] = j*dt

#plot das posicoes em funcao do tempo
for i in range(N-1):
    plt.plot(tempo, u[i])

plt.xlim(0, T)
#print(u[2])

#plt.show()
"""

tf = time.time()

print('Tempo de processamento (minutos): ', (tf-t0)/60)