# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import BEMmodel as BEM

""" ---- User defined variables: Wind Turbine ---- """
N = 100 # annulus spacing
mode = 'constant' # 2 possible modes: constant and cosinus distribution of annulus. 
                  # Write 'cosinus' or 'constant' if you want to change the distribution.
plotter = True # If you don't want to plot all the default plots please enter False

TSR = [6, 8, 10]

[results, CT, CP, Uinf, Radius] = BEM.ExecuteBEM(N, plotter, mode, TSR[0])
alpha = np.zeros([len(results[:,0]),3])
a = np.zeros([len(results[:,0]),3])
ap = np.zeros([len(results[:,0]),3])
phi = np.zeros([len(results[:,0]),3])
Ct = np.zeros([len(results[:,0]),3])
Cn = np.zeros([len(results[:,0]),3])
Cq = np.zeros([len(results[:,0]),3])
Circ = np.zeros([len(results[:,0]),3])

for i in range(0,3):
    [results, CT, CP, Uinf, Radius] = BEM.ExecuteBEM(N, plotter, mode, TSR[i])
    alpha[:,i] = results[:,7]
    phi[:,i] = results[:,6]
    a[:,i] = results[:,0]
    ap[:,i] = results[:,1]
    Ct[:,i] = results[:,4]
    Cn[:,i] = results[:,3]
    Cq[:,i] = results[:,8] 
    Circ[:,i] = results[:,5]   

""" Plots Section II.C """

# axial induction vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], alpha[:,0], 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], alpha[:,1], 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], alpha[:,2], 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\alpha$ [deg]', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Alpha_WT_Res.pdf', format='pdf', dpi=1000)

# phi vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], phi[:,0], 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], phi[:,1], 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], phi[:,2], 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\phi$ [deg]', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('phi_WT_Res.pdf', format='pdf', dpi=1000)


# axial induction a vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], a[:,0], 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], a[:,1], 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], a[:,2], 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$a$', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('a_WT_Res.pdf', format='pdf', dpi=1000)

# axial induction a' vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], ap[:,0], 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], ap[:,1], 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], ap[:,2], 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$a^,$', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('ap_WT_Res.pdf', format='pdf', dpi=1000)

# Cn vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], Cn[:,0]/(0.5*Uinf**2*Radius), 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], Cn[:,1]/(0.5*Uinf**2*Radius), 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], Cn[:,2]/(0.5*Uinf**2*Radius), 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$C_n$', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Cn_WT_Res.pdf', format='pdf', dpi=1000)

# Ct vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], Ct[:,0]/(0.5*Uinf**2*Radius), 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], Ct[:,1]/(0.5*Uinf**2*Radius), 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], Ct[:,2]/(0.5*Uinf**2*Radius), 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$C_t$', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Ct_WT_Res.pdf', format='pdf', dpi=1000)

# Cn vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], Cq[:,0]/(0.5*Uinf**2*Radius**2), 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], Cq[:,1]/(0.5*Uinf**2*Radius**2), 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], Cq[:,2]/(0.5*Uinf**2*Radius**2), 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$C_q$', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Cq_WT_Res.pdf', format='pdf', dpi=1000)

# Circulation vs. r/R
NBlades = 3
fig1 = plt.figure(figsize=(12, 6))
plt.plot(results[:,2], Circ[:,0]/(np.pi*Uinf**2/(NBlades*Uinf*TSR[0]/Radius)), 'r-', label=r'$\lambda=6$')
plt.plot(results[:,2], Circ[:,1]/(np.pi*Uinf**2/(NBlades*Uinf*TSR[1]/Radius)), 'b-', label=r'$\lambda=8$')
plt.plot(results[:,2], Circ[:,2]/(np.pi*Uinf**2/(NBlades*Uinf*TSR[2]/Radius)), 'g-', label=r'$\lambda=10$')
plt.grid()
plt.tick_params(axis='x', which='major', labelsize=20)
plt.tick_params(axis='y', which='major', labelsize=20)
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$\Gamma$', fontsize=20)
plt.legend(fontsize=20)
plt.savefig('Gamma_WT_Res.pdf', format='pdf', dpi=1000)
