# import necessary libraries
#from platypus import NSGAII, Problem, Real
import matplotlib.pyplot as plt
import numpy as np
import BEM_Opt as BEM_Opt
import BEMmodel as BEM

""" ---- Code variables ---- """
N = 100 # annulus spacing

""" ---- Optimization ---- """

# BACKGROUND: This optimization is based on the concept of optimizing the performance
# of the rotor for maximum Cp with two approaches:
# Propeller case: Ct = Free choice and Maximum Efficiency at all radial positions
    
CT0 = 0.12
CT_new = 0.12*(1+0.25**2) #smaller area
NB = 6
Radius = 0.7
Uinf = 60
Omega = 1200/60*2*np.pi
TSR = Omega*Radius/Uinf

Variables = [CT_new, NB, Radius, Uinf, TSR, Omega]
    
finish = False  

count=1
while (finish==False):
    
    [CP, CT, fnorm, ftan, gamma, chord, twist, r_R, phi, a, aline, vmag2] = BEM_Opt.ExecuteBEM_opt(N, Variables)
    count=count+1
    if (abs(CT0-CT)<0.00001):
        
        finish = True
        
    else:
        
        Variables[0]=Variables[0] + 0.05*(CT0-CT)
        
        if (count/100-np.int(count/100) == 0):
            display(CT-CT0)
        
"""" Baseline case """
mode = 'constant' # 2 possible modes: constant and cosinus distribution of annulus. 
                  # Write 'cosinus' or 'constant' if you want to change the distribution.
plotter = False # If you don't want to plot all the default plots please enter False

[results, CT, CP, Uinf, Radius] = BEM.ExecuteBEM(N, plotter, mode)
        
""" Plotter """

c_ba = (0.18-0.06*(r_R))*Radius # meters
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Chord distribution in Propeller', fontsize=20)
plt.plot(r_R,c_ba, 'k-', label=r'Baseline')
plt.plot(r_R, chord, 'r-', label=r'$Opt.$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$c(r/R)$ [m]', fontsize=20)
plt.savefig('Chord_distrib_PR.pdf', format='pdf', dpi=1000)

# Twist distribution
# Chord distribution
col_pitch = 46
t_ba = (-50*(r_R)+35+col_pitch) # degrees
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Twist distribution in Propeller', fontsize=20)
plt.plot(r_R,t_ba, 'k-', label=r'Baseline')
plt.plot(r_R, twist, 'r-', label=r'Opt.')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\theta(r/R)$ [deg]', fontsize=20)
plt.savefig('Twist_distrib_PR.pdf', format='pdf', dpi=1000)

# Fnorm distribution
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Normal force distribution in Propeller', fontsize=20)
plt.plot(results[:,2], results[:,3]/(0.5*Uinf**2*Radius), 'k-', label=r'$Baseline$')
plt.plot(r_R, fnorm/(0.5*Uinf**2*Radius), 'r-', label=r'Opt.')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$f_n(r/R)$ [N]', fontsize=20)
plt.savefig('fnorm_distrib_PR.pdf', format='pdf', dpi=1000)

# Circulation distribution
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Circulation distribution in Propeller', fontsize=20)
plt.plot(results[:,2], results[:,5]/(np.pi*Uinf**2/(NB*Uinf*TSR/Radius)), 'k-', label=r'Baseline')
plt.plot(r_R, gamma/(np.pi*Uinf**2/(NB*Uinf*TSR/Radius)), 'r-', label=r'Opt.')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\Gamma(r/R)$ [1/s]', fontsize=20)
plt.savefig('Gamma_PR.pdf', format='pdf', dpi=1000)



    

