import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import pandas as pd
import os
############     Dimensiones     ############
a = .1 #longitud
b = .158 #altura
e = 0.001 # espesor
############     Propiedades del material     ############
alfa = 1E-5 #dif. ter.
k = 30 #cond. term.
############     Parametros del Sistema     ############
T0 = 20 #temp. ini.
T_aire = 25 #temp. aire
v = 0.01167 #velocidad
L = .1 #distancia recorrida
h = 0 #coef. conv.
tmax = 100 #tiempo máximo
############     Fuente calor     ############
P = .656*150*18.2 #potencia
ahf = 0.00243
ahb = 0.00696
bh = 0.00296
rf = ahf/(ahf+ahb)
rb = ahb/(ahf+ahb)
############     Cálculo del modelo     ############
deltax = .001
deltay = .001
deltat = 0.01

x = np.linspace(0,a,int(a*1000+1))
y = np.linspace(0,b,int(b*1000+1))
X, Y = np.meshgrid(x,y)
t = np.linspace(0,tmax,int(tmax*100+1))
T = np.zeros((len(t),len(y),len(x)))
T[0,:,:] = T0

lamb = deltat*alfa/deltax**2
mu = deltat*alfa/deltay**2
print(lamb,mu)
############     Outputs     ############

print('{:-<45}\n{:^8}{:^5.2f}|{:^8}{:^5.2f}|\n{:-<45}'
        .format('','lambda:',lamb,'mu:',mu,''))

X0 = 0.05
Y0 = 0.02

sig = np.zeros(len(t))
sig[np.where(t<round(L/v,2))]=1
pos = np.zeros(len(t))
pos = Y0+v*t*sig
pos[int(round(L/v,2)*100+1)-1:] = pos[int(round(L/v,2)*100+1)-2]

def calc_2D_flux():
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
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
    # ax1.scatter(X,Y,c='k',s=.1)
    cmap = plt.get_cmap('jet', 200)
    cmap.set_over('grey')
    im = ax1.imshow([[]],cmap=cmap,interpolation='bilinear', extent=[0,a,0,b],origin='lower')

    cb = fig.colorbar(im, orientation='vertical',extend='max',ax=ax1)
    cb.ax.tick_params(labelsize=8)
    ax1.set_xlabel('Width (mm)')
    ax1.set_ylabel('Height (mm)', rotation=90)
    cb.set_label('Temp. ($^\circ$C)', fontsize=8)
    t_title = ax1.set_title('', fontsize=15)
    w, = ax1.plot([],[],'ko',mfc='none',markersize=10,markeredgewidth=2)
    w_ln, = ax1.plot([],[],'y',lw=10, solid_capstyle='round',alpha=.2)

    ax2.set_title('animation')
    ax2.plot(t,sig)
    ax2.set_ylim([-.5,1.5])
    l = ax2.axvline(x=0,c='r')

    def init():
        cb.set_ticks(np.linspace(20,1800,5))
        im.set_data(np.zeros((len(x),len(y))))
        w.set_data(X0,Y0)
        w_ln.set_data(X0,Y0)
        return im,w,w_ln

    def anim(i):
        im.set_clim(20,1800)
        im.set_data(T[i,:,:])
        w.set_data(X0,pos[i])
        w_ln.set_data([X0,X0],[Y0,pos[i]])
        l.set_xdata([t[i],t[i]])

        t_title.set_text('Time: {:4.1f} sec.'.format(t[i]))
        # if t[i]==0 or t[i]==5 or t[i]==10 or t[i]==15 or t[i]==20:
            # plt.savefig('welding_sim_{}_s.png'.format(t[i]))
        return im,w,w_ln

    if plotRealTime:
        an = FuncAnimation(fig, anim, frames=len(t), init_func=init,interval=1, repeat=False)
    else:
        anim(len(t)-1)
    plt.show()

def anim_2D_2(T,plotRealTime):
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
    # ax1.scatter(X,Y,c='k',s=.1)
    cmap = plt.get_cmap('jet', 200)
    cmap.set_over('grey')
    im = ax1.imshow([[]],cmap=cmap,interpolation='bilinear', extent=[0,a,0,b],origin='lower')

    cb = fig.colorbar(im, orientation='vertical',extend='max',ax=ax1)
    cb.ax.tick_params(labelsize=8)
    ax1.set_xlabel('Width (mm)')
    ax1.set_ylabel('Height (mm)', rotation=90)
    cb.set_label('Temp. ($^\circ$C)', fontsize=8)
    t_title = ax1.set_title('', fontsize=15)
    w, = ax1.plot([],[],'ko',mfc='none',markersize=10,markeredgewidth=2)
    w_ln, = ax1.plot([],[],'y',lw=10, solid_capstyle='round',alpha=.2)
    ax1.scatter([x[40],x[60],x[40],x[60]],[y[40],y[40],y[100],y[100]],c=['b','r', 'g', 'm'],s=100,edgecolors='black',linewidth=1.5)

    ax2.set_title('Control Point Temp.($^\circ$C)')
    l1, = ax2.plot([],[],'b')
    l2, = ax2.plot([],[],'g')
    ax2.axis([0, tmax+1, 0, 1000])

    def init():
        cb.set_ticks(np.linspace(20,1800,5))
        im.set_data(np.zeros((len(x),len(y))))
        w.set_data(X0,Y0)0
        w_ln.set_data(X0,Y0)
        return im,w,w_ln

    def anim(i):
        im.set_clim(20,1800)
        im.set_data(T[i,:,:])
        w.set_data(X0,pos[i])
        w_ln.set_data([X0,X0],[Y0,pos[i]])
        l1.set_data(t[:i],T_[:i,40,40])
        l2.set_data(t[:i],T_[:i,100,40])

        t_title.set_text('Time: {:4.1f} sec.'.format(t[i]))
        # if t[i]==0 or t[i]==5 or t[i]==10 or t[i]==15 or t[i]==20:
            # plt.savefig('welding_sim_{}_s.png'.format(t[i]))
        return im,w,w_ln

    if plotRealTime:
        an = FuncAnimation(fig, anim, frames=len(t), init_func=init,interval=1, repeat=False)
    else:
        anim(len(t)-1)
    plt.show()


def guardar_res(name_csv):
    datos_para_guardar = pd.DataFrame({'Tiempo':t,'Punto_1':T_[:,40,40],'Punto_2':T_[:,100,40]})
    datos_para_guardar.to_csv(os.path.join('data', '{}.csv'.format(name_csv)))

T_ = calc_2D_flux()
# anim_2D(T_,False)
anim_2D_2(T_,False)
####################################################################################
# guardar_res('resultados_sin_con')


