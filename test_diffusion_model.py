# This script simulates the passage of a charged particle through a positive
# ion plasma

## To retrieve the error plot from the excellent article by Delforge & co. the
## following values were used, sorted by decreasing interparticle distance to 
## Debye length ratios. But beware : the lower the ratio the longer the computation.
## The last sets of parameters take several hours to compute, the first ones are way faster
    
    
##  n_0 = 1e22 ; T = k_B/e*(1e3) ; N = int(1e6) ; dt = 1e-17 ( MHD Generator plasma )

##  n_0 = 1e25 ; T = k_B/e*(1e5) ; N = int(1e5) ; dt = 1e-17 ( Laser-produced plasma )

##  n_0 = 1e18 ; T = k_B/e*(1e3) ; N = int(1e5) ; dt = 1e-17 ( Alcaline plasma )

##  n_0 = 1e26 ; T = k_B/e*(1e6) ; N = int(1e5) ; dt = 1e-17 ( Thermonuclear plasma )

##  n_0 = 1e21 ; T = k_B/e*(1e5) ; N = int(1e5) ; dt = 1e-16 ( High intensity discharge plasma )

##  n_0 = 1e17 ; T = k_B/e*(1e4) ; N = int(1e5) ; dt = 1e-15 ( Low intensity discharge plasma )

##  n_0 = 3e11 ; T = k_B/e*(1e3) ; N = int(1e5) ; dt = 1e-14 ( High ionosphere plasma )


from classes import *
import numpy as np
from graphes import *
from champs import E_plasma as f_E



### Numerical constants ###


Nt = 1000000 ## number of time steps

dt = 1e-17 ## time step (s)

DT = dt/2


### Physics constants ###


e = 1.602e-19 ## elementary charge (C)

u = 1.66e-27 ## atomic mass unit (kg)

epsilon = 8.85e-12 ## vacuum dielectric permittivity (A2.s4/m3/kg)

k_B = 1.380649e-23 ## Boltzmann constant (m2.kg./K/s2)

c = 299792458 ## velocity of light in vacuum (m/s)



q = -1*e ## incident particle charge (C) (-1 for e-, +2 for alpha)

m = (5.49e-4)*u ## incident particle mass (kg) (5.49e-4 for e-, 4.00153 for alpha)

Z = +1 ## heavy nucleus atomic number Z



T = k_B*(1e3)/e ## ion plasma temperature (eV)

n_0 = 1e22 ## ion plasma density (m-3)

d_0 = n_0**(-1/3) ## ion plasma interparticle distance (m)

debye = np.sqrt(epsilon*T/(n_0*e*Z**2)) ## Debye length (m)

print('\n Ratio of interparticle distance by Debye length = ' + str(d_0/debye))



###############################################################################
#                           INITIALIZATION                                    #
###############################################################################


# position vector (m)

X = Vector(-1*debye, 8*d_0/10, 7*d_0/10)


# velocity vector (m/s)

V = Vector((1/100)*c, 0, 0)

# relativistic impulsion vector (kg.m/s)

P = U(V).produit(m)


# Trajectory parameters storing arrays

traj_x, traj_y, traj_z = [], [], []

V_x, V_y, V_z = [], [], []

P_x, P_y, P_z = [], [], []

Ec, temps = [], []


# Fields

B = Vector(0, 0, 0)

E = Vector(0, 0, 0)

E0 = 1e6 ## stationnary electric field for the Drude model test

# Instants of electron-ion interaction variables

I = 0

instants = []


###############################################################################
#                          PARTICLE PROPAGATION                               #
###############################################################################



for p in range(Nt):
    
    Xinter = X.plus(V.produit(DT))
    
    E = f_E(Z, Xinter, d_0, debye) ## storing array for the electric field and the interaction counter
    phase = step(Xinter, V, B, E, q, m, dt) # storing array for [X,V]
    
    X = phase[0] # position
    
    V = phase[1] # velocity
    
    P = U(V).produit(m) # impulsion
    
    
    traj_x.append(X.x*100)
    traj_y.append(X.y*100)
    traj_z.append(X.z*100)
    
    V_x.append(V.x)
    V_y.append(V.y)
    V_z.append(V.z)
    
    P_x.append(P.x)
    P_y.append(P.y)
    P_z.append(P.z)
    
    
    Ec.append( m*produit_scalaire(V, V)/(2*e )) # kinetic energy
    
    temps.append(dt*p*10**12) # time



###############################################################################
#                                    PLOTS                                    #
###############################################################################



traj_x = np.array(traj_x)/(100*d_0)
MAX_x = int(max(traj_x)//1)
MIN_x = int(min(traj_x)//1)

traj_y = np.array(traj_y)/(100*d_0)
MAX_y = int(max(traj_y)//1)
MIN_y = int(min(traj_y)//1)

traj_z = np.array(traj_z)/(100*d_0)
MAX_z = int(max(traj_z)//1)
MIN_z = int(min(traj_z)//1)



plot_2d_trajs(traj_x, traj_y, "(x,y) Trajectory", r"$\frac{x}{d_0}$", r"$\frac{y}{d_0}$")

plt.xlim([-2, MAX_x+1])
plt.ylim([MIN_y-1, MAX_y+1])

# listex = []
# listey = []

# for i in range(MAX_x + 2):
#     for j in range(MIN_y,MAX_y+2):
#         listex.append(i)
#         listey.append(j)

# plt.plot(np.array(listex), np.array(listey), '.', color = 'Black', label = 'Nucleus')
# plt.legend()

plt.show()



# plot_2d_trajs(traj_z, traj_y, "(z,y) Trajectory", r"$\frac{z}{d_0}$", r"$\frac{y}{d_0}$")

# plt.xlim([MIN_z-1, MAX_z+1])
# plt.ylim([MIN_y-1, MAX_y+1])

# listex = []
# listey = []

# for i in range(MIN_z,MAX_z+2):
#     for j in range(MIN_y,MAX_y+2):
#         listex.append(i)
#         listey.append(j)

# plt.plot(np.array(listex), np.array(listey), '.', color = 'Black', label = 'Nucleus')
# plt.legend()

# plt.show()



# plot_2d_trajs(traj_x, traj_z, "(x,z) Trajectory", r"$\frac{x}{d_0}$", r"$\frac{z}{d_0}$")

# plt.xlim([-2, MAX_x+1])
# plt.ylim([MIN_z-1, MAX_z+1])

# listex = []
# listey = []

# for i in range(MAX_x+2):
#     for j in range(MIN_z,MAX_z+2):
#         listex.append(i)
#         listey.append(j)

# plt.plot(np.array(listex), np.array(listey), '.', color = 'Black', label = 'Nucleus')
# plt.legend()

# plt.show()

plot_3d(traj_x, traj_y, traj_z, "3D Trajectory", r"$\frac{x}{d_0}$", r"$\frac{y}{d_0}$", r"$\frac{z}{d_0}$")

plot_2d_graphs(temps, Ec, "Kinetic Energy", "t [ps]", "Ec [eV]")

E_moy = sum(Ec)/len(Ec)

plt.plot(temps, np.array([E_moy for k in range(len(temps))]), color = 'Red', label = 'Average value')
plt.legend()
plt.show()

# plot_2d_graphs(temps, P_x, "Longitudinal impulsion", "t [ps]", r"$p_{||}$ [kg.m.s$^{-1}$]")
# plt.show()


tps = []
dtP = []
P2 = []

for i in range(len(temps)-1):
        
    tps.append(temps[i])
    dtP.append((P_x[i+1]-P_x[i])/dt)
    P2.append(2*dtP[i]/(P_x[i+1]+P_x[i]))


# plot_2d_graphs(tps, dtP, "Longitudinal impulsion derivative", "t [ps]", r"$\frac{dp_{||}}{dt}$ [kg.m.s$^{-2}$]")
# plt.yscale('log')
# plt.show()

f_moy = abs(sum(P2)/len(P2)) # average collision frequency (Hz)

plot_2d_graphs(tps, P2, "Collision frequency", "t [ps]", r"$\frac{1}{p_{||}}\frac{dp_{||}}{dt}$ [s$^{-1}$]")
plt.plot(tps, np.array([f_moy for i in range(len(tps))]), color = 'Red', label = 'Average value')
plt.yscale('log')
plt.legend()
plt.show()

bc = abs(q)*Z*e/(4*np.pi*epsilon*m*V_x[0]**2)

f_coulomb = n_0*V_x[0]**2*4*np.pi*bc**2*np.log((1+(bc/debye)**2)**0.5/(bc/debye))/(abs(sum(V_x)/len(V_x))) # Coulomb collision frequency (Hz)

print('\n Relative error on the Coulomb collision frequency = ' + str(abs(f_moy-f_coulomb)/f_coulomb))

fE = abs(E_moy-Ec[0])/(Ec[0]*Nt*dt)

Vx_moy = sum(V_x)/len(V_x)
fC = abs(Vx_moy-V_x[0])/(V_x[0]*Nt*dt)

print('\n ratio of kinetic by longitudinal collision frequency = ' + str(fE/fC))

