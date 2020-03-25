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

WTOpt = True # Select what you want to optimize
PROpt = False

if WTOpt == True:
    
    CT = 0.75*(1+0.2**2) #Corrected value, root at r=0.2, Area is lower
    NB = 3
    Radius = 50
    Uinf = 10
    TSR = 8
    Omega = Uinf*TSR/Radius
    
    Variables = [CT, NB, Radius, Uinf, TSR, Omega]
     
if PROpt == True:
    
    NB = 6
    Radius = 0.7
    Uinf = 60
    Omega = 1200/60*2*np.pi
    TSR = Omega*Radius/Uinf
    
    Variables = [NB, Radius, Uinf, TSR, Omega]
    
finish = False  
CT0 = 0.75 

while (finish==False):
    
    [CP, CT, fnorm, ftan, gamma, chord, twist, r_R, phi, a, aline, vmag2] = BEM_Opt.ExecuteBEM_opt(PROpt, WTOpt, N, Variables)

    if (abs(CT0-CT)<0.001):
        
        finish = True
        
    else:
        
        Variables[0]=Variables[0] + 0.05*(CT0-CT)
        display(CT-CT0)
        
""" Plotter """

# Chord distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Chord distribution in Propeller', fontsize=20)
plt.plot(r_R, chord, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$c(r/R)$ [m]', fontsize=20)
plt.savefig('Chord_distrib_WT.pdf', format='pdf', dpi=1000)

# Twist distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Twist distribution in Propeller', fontsize=20)
plt.plot(r_R, twist, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\theta(r/R)$ [deg]', fontsize=20)
plt.savefig('Twist_distrib_WT.pdf', format='pdf', dpi=1000)

# Fnorm distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Normal force distribution in Propeller', fontsize=20)
plt.plot(r_R, fnorm, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$f_n(r/R)$ [N]', fontsize=20)
plt.savefig('fnorm_distrib_WT.pdf', format='pdf', dpi=1000)

# Circulation distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Circulation distribution in Propeller', fontsize=20)
plt.plot(r_R, gamma, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\Gamma(r/R)$ [1/s]', fontsize=20)
plt.savefig('Gamma_WT.pdf', format='pdf', dpi=1000)

# Inflow angle distribution
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Inflow distribution in Propeller', fontsize=20)
plt.plot(r_R, a, 'r-')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$a(r/R)$ [1/s]', fontsize=20)
plt.savefig('a_WT.pdf', format='pdf', dpi=1000)

    

