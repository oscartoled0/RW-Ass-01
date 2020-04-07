# import necessary libraries
#from platypus import NSGAII, Problem, Real
import matplotlib.pyplot as plt
import numpy as np
import BEM_Opt as BEM_Opt
import math as m

""" ---- Code variables ---- """
N = 100 # annulus spacing

""" ---- Optimization ---- """

# BACKGROUND: This optimization is based on the concept of optimizing the performance
# of the rotor for maximum Cp with two approaches:
# Propeller case: Ct = Free choice and Maximum Efficiency at all radial positions
    
CT0 = 0.12
CT_new = 0.12*(1+0.25**2) #smaller area
Radius = 0.7
NB = 6
Uinf = 60
Omega = 1200/60*2*np.pi
TSR = Omega*Radius/Uinf
pinf = 94200
rho = 1

Variables = [CT_new, NB, Radius, Uinf, TSR, Omega]
    
finish = False  

count=1
while (finish==False):
    
    [CP, CT, fnorm, ftan, gamma, chord, twist, r_R, phi, a, aline, vmag2, r_ext] = BEM_Opt.ExecuteBEM_opt(N, Variables)
    count=count+1
    if (abs(CT0-CT)<0.00001):
        
        finish = True
        
    else:
        
        Variables[0]=Variables[0] + 0.05*(CT0-CT)
        
        if (count/100-np.int(count/100) == 0):
            display(CT-CT0)
        
""" Plotter """

# Chord distribution
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Chord distribution in Propeller', fontsize=20)
plt.plot(r_R, chord, 'r-')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$c(r/R)$ [m]', fontsize=20)
plt.savefig('Chord_distrib_PR.pdf', format='pdf', dpi=1000)

# Twist distribution
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Twist distribution in Propeller', fontsize=20)
plt.plot(r_R, twist, 'r-')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\theta(r/R)$ [deg]', fontsize=20)
plt.savefig('Twist_distrib_PR.pdf', format='pdf', dpi=1000)

# Fnorm distribution
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Normal force distribution in Propeller', fontsize=20)
plt.plot(r_R, fnorm, 'r-')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$f_n(r/R)$ [N]', fontsize=20)
plt.savefig('fnorm_distrib_PR.pdf', format='pdf', dpi=1000)

# Circulation distribution
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Circulation distribution in Propeller', fontsize=20)
plt.plot(r_R, gamma, 'r-')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\Gamma(r/R)$ [1/s]', fontsize=20)
plt.savefig('Gamma_PR.pdf', format='pdf', dpi=1000)

# Inflow angle distribution
fig1 = plt.figure(figsize=(12, 6))
# plt.title(r'Inflow distribution in Propeller', fontsize=20)
plt.plot(r_R, a, 'r-')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$a(r/R)$ [1/s]', fontsize=20)
plt.savefig('a_PR.pdf', format='pdf', dpi=1000)

Uw = np.zeros(a[:].size)
Ud = np.zeros(a[:].size)

for i in range(a.size):
    Uw[i] = m.sqrt((Uinf*(1+2*a[i]))**2 + (Omega*r_R[i]*Radius*aline[i])**2)
    Ud[i] = m.sqrt((Uinf*(1+a[i]))**2 + (Omega*r_R[i]*Radius*aline[i])**2)
    
Uw2 = Uinf*(1+2*a[:])
Ud2 = Uinf*(1+a[:])
    
Hinf = np.zeros(len(Uw))
H1mas = np.zeros(len(Uw))
H1menos = np.zeros(len(Uw))
Hminf = np.zeros(len(Uw))
for i in range(len(Uw)):
    Hinf[i] = pinf/rho + 0.5*Uinf**2
    H1mas[i] = pinf/rho + 0.5*Uinf**2
    H1menos[i] = pinf/rho + 0.5*Uw[i]**2
    Hminf[i] = pinf/rho + 0.5*Uw[i]**2

r_inft = np.zeros(len(r_ext))    
r_minft = np.zeros(len(r_ext))  
r_Rinft = np.zeros(len(r_R))    
r_Rminft = np.zeros(len(r_R))  
r_inft[0] = r_ext[0]  
r_minft[0] = r_ext[0]  
for i in range(1,len(r_ext)):
    r_inft[i] = m.sqrt(((r_ext[i]**2-r_ext[i-1]**2)*Ud2[i-1]/(Uinf)) + r_inft[i-1]**2)
    r_minft[i] = m.sqrt(((r_ext[i]**2-r_ext[i-1]**2)*Ud2[i-1]/(Uw2[i-1])) + r_minft[i-1]**2)
    r_Rinft[i-1] = (r_inft[i] + r_inft[i-1])/2
    r_Rminft[i-1] = (r_minft[i] + r_minft[i-1])/2


fig1 = plt.figure(figsize=(18, 8))
plt.plot((Hinf-500)/1000, r_Rinft, 'c-', label = r'$x/x_d=\infty$')
plt.plot(H1mas/1000, r_R, 'b-', label = r'$x/x_d=1^+$')
plt.plot(H1menos/1000, r_R, 'g-', label = r'$x/x_d=1^-$')
plt.plot((Hminf+500)/1000, r_Rminft, 'r-', label = r'$x/x_d=-\infty$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel(r'$H$ [kPa]', fontsize=20)
plt.ylabel(r'r/R', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('H1menosprop.pdf', format='pdf', dpi=1000)
    

