import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import numpy as np
# import pandas as pd
# import os
############################################################
#                             INPUTS
############################################################
############     Dimensiones     ############
a = .1 #longitud
b = .158 #altura
e = 0.001 # espesor
X0 = 0.05 #posición inicial x
Y0 = 0.02 #posición inicial y
L = .1 #distancia recorrida
############     Parametros del Sistema     ############
v = 0.01167 #velocidad
tmax = 50 #tiempo máximo
T0 = 20 #temp. ini.
T_aire = 25 #temp. aire
############     Propiedades del material     ############
k = 30 #cond. term.
alfa = 1E-5 #dif. ter.
h = 0 #coef. conv.
############################################################
#                             Cálculos
############################################################
deltax = .001
deltay = .001
deltat = 0.01
############################################################
x = np.linspace(0,a,int(a*1000+1))
y = np.linspace(0,b,int(b*1000+1))
X, Y = np.meshgrid(x,y)
t = np.linspace(0,tmax,int(tmax*100+1))
T = np.zeros((len(t),len(y),len(x)))
T[0,:,:] = T0
############     Fuente calor     ############
P = .656*150*18.2 #potencia
ahf = 0.00243
ahb = 0.00696
bh = 0.00296
rf = ahf/(ahf+ahb)
rb = ahb/(ahf+ahb)
############################################################
lamb = deltat*alfa/deltax**2
mu = deltat*alfa/deltay**2
############################################################
sig = np.zeros(len(t))
sig[np.where(t<round(L/v,2))]=1
pos = np.zeros(len(t))
pos = Y0+v*t*sig
pos[int(round(L/v,2)*100+1)-1:] = pos[int(round(L/v,2)*100+1)-2]

def calculo():
    for j in range(len(t)-1):
        Z =(6*rf*P)/(ahf*bh*np.pi)*np.exp(-3*((X-X0)**2/(ahf**2)+(Y-pos[j])**2/(bh**2)))*sig[j]
# centro
        T[j+1,1:-1,1:-1] = ((1-2*lamb-2*mu)*T[j,1:-1,1:-1]
                            +lamb*(T[j,1:-1,2:]+T[j,1:-1,:-2])
                            +mu*(T[j,2:,1:-1]+T[j,:-2,1:-1])
                            +Z[1:-1,1:-1]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,1:-1,1:-1]))
#linea vertical izq.
        T[j+1,1:-1,0] = ((1-2*lamb-2*mu)*T[j,1:-1,0]
                        +2*lamb*T[j,1:-1,1]
                        +mu*(T[j,2:,0]+T[j,:-2,0])
                        +Z[1:-1,0]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,1:-1,0]))
#linea vertical der.
        T[j+1,1:-1,-1] = ((1-2*lamb-2*mu)*T[j,1:-1,-1]
                         +2*lamb*T[j,1:-1,-2]
                         +mu*(T[j,2:,-1]+T[j,:-2,-1])
                         +Z[1:-1,-1]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,1:-1,-1]))
#linea horizontal de abajo
        T[j+1,0,1:-1] =  ((1-2*lamb-2*mu)*T[j,0,1:-1]
                         +lamb*(T[j,0,2:]+T[j,0,:-2])
                         +2*mu*T[j,1,1:-1]
                         +Z[0,1:-1]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,0,1:-1]))
#linea horizontal de arriba
        T[j+1,-1,1:-1] = ((1-2*lamb-2*mu)*T[j,-1,1:-1]
                         +lamb*(T[j,-1,2:]+T[j,-1,:-2])
                         +2*mu*T[j,-2,1:-1]
                         +Z[-1,1:-1]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,-1,1:-1]))
# UR_C (upper/right)
        T[j+1,-1,-1] = ((1-2*lamb-2*mu)*T[j,-1,-1]
                         +2*lamb*T[j,-1,-2]
                         +2*mu*T[j,-2,-1]
                         +Z[-1,-1]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,-1,-1]))
# UL_R upper/left
        T[j+1,-1,0] = ((1-2*lamb-2*mu)*T[j,-1,0]
                         +2*lamb*T[j,-1,1]
                         +2*mu*T[j,-2,0]
                         +Z[-1,0]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,-1,0]))
# BR_C (bottom/right)
        T[j+1,0,-1] = ((1-2*lamb-2*mu)*T[j,0,-1]
                         +2*lamb*T[j,0,-2]
                         +2*mu*T[j,1,-1]
                         +Z[0,-1]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,0,-1]))
# BL_C (bottom/left)
        T[j+1,0,0] = ((1-2*lamb-2*mu)*T[j,0,0]
                         +2*lamb*T[j,0,1]
                         +2*mu*T[j,1,0]
                         +Z[0,0]*deltat*alfa/(k*e)+(2*h*deltat*alfa/(k*e))*(T_aire-T[j,0,0]))
    return T


def anim_2D(T,plotRealTime):
    fig, ax = plt.subplots(figsize=(8,8))
    cmap = plt.get_cmap('jet', 200)
    cmap.set_over('grey')
    im = ax.imshow([[]],cmap=cmap,interpolation='bilinear', extent=[0,a,0,b],origin='lower')

    cb = fig.colorbar(im, orientation='vertical',extend='max',ax=ax)
    cb.ax.tick_params(labelsize=8)
    ax.set_xlabel('Width (mm)',fontsize=12)
    ax.set_ylabel('Height (mm)',rotation=90,fontsize=12)
    cb.set_label('Temp. ($^\circ$C)', fontsize=12)
    t_title = ax.set_title('', fontsize=15)
    w, = ax.plot([],[],'ko',mfc='none',markersize=10,markeredgewidth=2)

    def init():
        cb.set_ticks(np.linspace(20,1800,5))
        im.set_data(np.zeros((len(x),len(y))))
        w.set_data(X0,Y0)
        return im,w

    def anim(i):
        im.set_clim(20,1800)
        im.set_data(T[i,:,:])
        w.set_data(X0,pos[i])
        t_title.set_text('Time: {:4.1f} sec.'.format(t[i]))

        return im,w

    if plotRealTime:
        an = FuncAnimation(fig, anim, frames=len(t), init_func=init,interval=1, repeat=False)
    else:
        anim(len(t)-1)
    plt.show()


def points_temp(T,XcontrolP,YcontrolP,plotRealTime):
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
    colors = cm.jet(np.linspace(0, 1, len(XcontrolP)))
    ax1.scatter(X,Y,c='k',s=.01)
    ax1.set_xlabel('Width (mm)')
    ax1.set_ylabel('Height (mm)', rotation=90)
    ax1.axvline(x=a/2,color='k',linewidth=.8,linestyle='--')
    ax1.scatter(XcontrolP,YcontrolP,c=colors,s=100,edgecolors='black',linewidth=.5)

    ax2.set_title('Control Point Temp.($^\circ$C)')
    ax1.axis([0, a, 0, b])
    lines = [ax2.plot([0,1],[0,YcontrolP[l]],c = colors[l]) for l in range(len(XcontrolP))]

    ax2.axis([0, tmax+1, 0, 1800])

    def init():
        [l[0].set_data(0,0) for j,l in enumerate(lines)]

    def anim(i):
        [l[0].set_data(t[:i],T_[:i,int(YcontrolP[j]*1000),int(XcontrolP[j]*1000)]) for j,l in enumerate(lines)]

    if plotRealTime:
        an = FuncAnimation(fig, anim, frames=len(t), init_func=init,interval=1, repeat=False)
    else:
        anim(len(t)-1)
    plt.show()

T_ = calculo()
# anim_2D(T_,True)
points_temp(T_,[0.043,0.045,0.06],[0.05,0.1,0.075],False)
