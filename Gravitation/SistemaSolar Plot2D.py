import numpy as np
import matplotlib.pyplot as plt

# Parámetros
N = 9
G = 8.88743163e-10
Nt = 250000
dt = 0.25

M = np.zeros((N))
M[0] = 332970.529136    # Sol.
M[1] = 0.0552735261     # Mercurio.
M[2] = 0.814997513      # Venus.
M[3] = 1                # La Tierra.
M[4] = 0.107446849      # Marte.
M[5] = 317.828133       # Júpiter.
M[6] = 95.16            # Saturno.
M[7] = 14.54            # Urano.
M[8] = 17.15            # Neptuno.

# Matrices de output.
X = np.zeros((N,Nt,3))
V = np.zeros((N,3))

# Condiciones iniciales.
X[0][0] = [0.0, 0.0, 0.0]
X[1][0] = [0.466700788, 0.0, 0.056916773]
X[2][0] = [0.728237503, 0.0, 0.043121576]
X[3][0] = [1.016713884, 0.0, -2.71676e-07]
X[4][0] = [1.666015896, 0.0, 0.053774992]
X[5][0] = [5.454635139, 0.0, 0.124169614]
X[6][0] = [10.05033838, 0.0, 0.435934741]
X[7][0] = [20.09599544, 0.0, 0.270987774]
X[8][0] = [30.32823783, 0.0, 0.936783848]

V[0] = [0.0, 0.0, 0.0]
V[1] = [0.0, 0.022443034, 0.0]
V[2] = [0.0, 0.020089909, 0.0]
V[3] = [0.0, 0.016917345, 0.0]
V[4] = [0.0, 0.012689974, 0.0]
V[5] = [0.0, 0.007185195, 0.0]
V[6] = [0.0, 0.005278104, 0.0]
V[7] = [0.0, 0.003745623, 0.0]
V[8] = [0.0, 0.003110241, 0.0]

# Evolución temporal.

# t - tiempo
# a - masa que se mueve
# b - otras masas
# j - número de k\q del RK4
# i - la componente del vector

def deltaX(t, a, b, j):
    W = X[b][t]-X[a][t]
    
    if j in range(1, 3):
        W = W + (1/2)*dt*K[j-1][b] - (1/2)*dt*K[j-1][a]
    elif j == 3:
        W = W + dt*K[j-1][b] - dt*K[j-1][a]
    
    return W
        
def dist(t, a, b, j):
    return ((deltaX(t, a, b, j)[0])**2+(deltaX(t, a, b, j)[1])**2+(deltaX(t, a, b, j)[2])**2)**(1/2)

def g(t, a, b, j):
    return G*(M[b]/(dist(t, a, b, j))**3)*(deltaX(t, a, b, j))

for t in range(Nt-1):
    Q = np.zeros((4,N,3))
    K = np.zeros((4,N,3))
    
    for j in range(4):
        for a in range(N):
            for b in range(N):
                if b == a:
                    continue
                else:
                    Q[j][a] = Q[j][a] + g(t, a, b, j)
            
            if j == 0:
                K[j][a] = V[a]
            elif j in range(1, 3):
                K[j][a] = V[a] + (dt/2)*Q[j][a]
            else:
                K[j][a] = V[a] + dt*Q[j][a]
    
    for a in range(N):
        V[a] = V[a] + (dt/6)*(Q[0][a]+2*Q[1][a]+2*Q[2][a]+Q[3][a])
        X[a][t+1] = X[a][t] + (dt/6)*(K[0][a]+2*K[1][a]+2*K[2][a]+K[3][a])
        
        if t%7500 == 0:
            print(t, V[a])
            print(t, X[a][t+1])
        
# Plot.
'''
for a in range(N):
    ax = plt.axes(projection='3d')
    ax.plot(X[a][0], X[a][1], X[a][2], color='black')
plt.show()
'''
Xplot = np.zeros((N,3,Nt))
for t in range(Nt):
    for a in range(N):
        for i in range(3):
            Xplot[a][i][t] = X[a][t][i]

plt.axes()
plt.figure(figsize=(10, 10))
plt.xlim((-35.0, 35.0))
plt.ylim((-35.0, 35.0))
plt.xlabel('x (UA)', fontsize = 14)
plt.ylabel('y (UA)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

s= r'$ Sol$'
plt.plot(Xplot[0][0], Xplot[0][1], color = 'k',label = s)

m= r'$ Mercurio$'
plt.plot(Xplot[1][0], Xplot[1][1],'#E50000', label=m)

v= r'$ Venus$'
plt.plot(Xplot[2][0], Xplot[2][1],'#370EFF', label=v)  

t= r'$ Tierra$'
plt.plot(Xplot[3][0], Xplot[3][1],'#895353', label=t) 

mm= r'$ Marte$'
plt.plot(Xplot[4][0], Xplot[4][1],'#10BD0B', label=mm)

j= r'$ Júpiter$'
plt.plot(Xplot[5][0], Xplot[5][1],'#FF8C00', label=j)

s= r'$ Saturno$'
plt.plot(Xplot[6][0], Xplot[6][1],'#DC143C', label=s)  

u= r'$ Urano$'
plt.plot(Xplot[7][0], Xplot[7][1],'#002BFF', label=u) 

n= r'$ Neptuno$'
plt.plot(Xplot[8][0], Xplot[8][1],'#00CDFF', label=n)
plt.legend(loc='upper right', prop={'size': 14})