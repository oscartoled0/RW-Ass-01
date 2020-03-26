# import necessary libraries
#from platypus import NSGAII, Problem, Real
import matplotlib.pyplot as plt
import numpy as np
import BEM_Opt as BEM_Opt

""" ---- Code variables ---- """
N = 100 # annulus spacing

""" ---- Optimization ---- """

# BACKGROUND: This optimization is based on the concept of optimizing the performance
# of the rotor for maximum Cp with two approaches:
# Wind Turbine case: Ct = 0.75 and Maximum Efficiency at all radial positions
# Propeller case: Ct = Free choice and Maximum Efficiency at all radial positions
    
CT_new = 0.75*(1+0.2**2) #Corrected value, root at r=0.2, Area is lower
NB = 3
Radius = 50
Uinf = 10
TSR = 8
Omega = Uinf*TSR/Radius
TSR = [6, 8, 10]
             
finish = False  
CT0 = 0.75 

Variables = [CT_new, NB, Radius, Uinf, TSR[0]]
[CP, CT, fnorm, ftan, gamma, chord, twist, r_R, phi, a, aline, vmag2] = BEM_Opt.ExecuteBEM_opt(N, Variables)


fnorm = np.zeros([len(r_R),3])
ftan = np.zeros([len(r_R),3])
gamma = np.zeros([len(r_R),3])
chord = np.zeros([len(r_R),3])
twist = np.zeros([len(r_R),3])
phi = np.zeros([len(r_R),3])
a = np.zeros([len(r_R),3])
aline = np.zeros([len(r_R),3])
vmag2 = np.zeros([len(r_R),3])
CP = np.zeros(3)
CT = np.zeros(3)
    
for i in range(0,3):
    
    finish = False
    while (finish==False):
        
        Variables = [CT_new, NB, Radius, Uinf, TSR[i]]
        [CP[i], CT[i], fnorm[:,i], ftan[:,i], gamma[:,i], chord[:,i], twist[:,i], r_R, phi[:,i], a[:,i], aline[:,i], vmag2[:,i]] = BEM_Opt.ExecuteBEM_opt(N, Variables)

        if (abs(CT0-CT[i])<0.0001):
            
            finish = True
            
        else:
            
            CT_new = CT_new + 0.05*(CT0-CT[i])
        
""" Plotter """

# Chord distribution
fig1 = plt.figure(figsize=(12, 6))
plt.plot(r_R, chord[:,0], 'r-', label = r'$\lambda=6$')
plt.plot(r_R, chord[:,1], 'b-', label = r'$\lambda=8$')
plt.plot(r_R, chord[:,2], 'g-', label = r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$c(r/R)$ [m]', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Chord_distrib_WT.pdf', format='pdf', dpi=1000)

# Twist distribution
fig1 = plt.figure(figsize=(12, 6))
plt.plot(r_R, twist[:,0], 'r-',  label = r'$\lambda=6$')
plt.plot(r_R, twist[:,1], 'b-', label = r'$\lambda=8$')
plt.plot(r_R, twist[:,2], 'g-', label = r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\theta(r/R)$ [deg]', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Twist_distrib_WT.pdf', format='pdf', dpi=1000)

# Fnorm distribution
fig1 = plt.figure(figsize=(12, 6))
plt.plot(r_R, fnorm[:,0]/(0.5*Uinf**2*Radius), 'r-',  label = r'$\lambda=6$')
plt.plot(r_R, fnorm[:,1]/(0.5*Uinf**2*Radius), 'b-', label = r'$\lambda=8$')
plt.plot(r_R, fnorm[:,2]/(0.5*Uinf**2*Radius), 'g-', label = r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$C_n(r/R)$ ', fontsize=20)
plt.savefig('fnorm_distrib_WT.pdf', format='pdf', dpi=1000)
plt.legend(fontsize=20)

# Circulation distribution
fig1 = plt.figure(figsize=(12, 6))
plt.plot(r_R, gamma[:,0]/(np.pi*Uinf**2/(NB*Uinf*TSR[0]/Radius)), 'r-', label = r'$\lambda=6$')
plt.plot(r_R, gamma[:,1]/(np.pi*Uinf**2/(NB*Uinf*TSR[0]/Radius)), 'b-', label = r'$\lambda=8$')
plt.plot(r_R, gamma[:,2]/(np.pi*Uinf**2/(NB*Uinf*TSR[0]/Radius)), 'g-', label = r'$\lambda=10$')
plt.xlabel('r/R', fontsize=20)
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.ylabel(r'$\Gamma(r/R)$', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Gamma_WT.pdf', format='pdf', dpi=1000)

# Inflow angle distribution
fig1 = plt.figure(figsize=(12, 6))
plt.plot(r_R, phi[:,0], 'r-', label = r'$\lambda=6$')
plt.plot(r_R, phi[:,1], 'b-', label = r'$\lambda=8$')
plt.plot(r_R, phi[:,2], 'g-', label = r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\phi(r/R)$ [deg]', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('phi_WT.pdf', format='pdf', dpi=1000)

# Inflow angle distribution
fig1 = plt.figure(figsize=(12, 6))
plt.plot(r_R, a[:,0], 'r-', label = r'$\lambda=6$')
plt.plot(r_R, a[:,1], 'b-', label = r'$\lambda=8$')
plt.plot(r_R, a[:,2], 'g-', label = r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$a(r/R)$ ', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('a_WT.pdf', format='pdf', dpi=1000)

    

