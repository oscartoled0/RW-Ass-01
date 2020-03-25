# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import BEMmodel as BEM

""" ---- User defined variables ---- """
PR = True
WT = False
N = 100 # annulus spacing
mode = 'constant' # 2 possible modes: constant and cosinus distribution of annulus. 
                  # Write 'cosinus' or 'constant' if you want to change the distribution.
plotter = False # If you don't want to plot all the default plots please enter False

[results, CT, CP, Uinf, Radius] = BEM.ExecuteBEM(PR, WT, N, plotter, mode)

""" Plots Section II.C """

# Alpha and inflow angle vs. r/R
fig, ax1 = plt.subplots(figsize=(12,6))
color = 'red'
ax1.set_xlabel('r/R', fontsize=20)
ax1.set_ylabel(r'$\alpha$ [deg]', color=color, fontsize=20)
ax1.plot(results[:,2], results[:,7], color=color, label=r'$\alpha$')
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid()
ax1.set_title('Angle of attack radial distribution', fontsize=20)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'blue'
ax2.set_ylabel(r'$\phi$ [deg]', color=color, fontsize=20)  # we already handled the x-label with ax1
ax2.plot(results[:,2], results[:,6], color=color, label=r'$\phi$')
ax2.tick_params(axis='y', labelcolor=color)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('Alpha_Inflow_lamb10.pdf', format='pdf', dpi=1000)

# axial induction vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.title('Axial and tangential induction radial distribution', fontsize=20)
plt.plot(results[:,2], results[:,0], 'r-', label=r'$a$')
plt.plot(results[:,2], results[:,1], 'b-', label=r'$a^,$')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$a$ and $a^,$', fontsize=20)
plt.legend(fontsize=20)
#plt.savefig('Axial_Tang_Induct_lamb10.pdf', format='pdf', dpi=1000)

# Cn Ct vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Normal and tagential force, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$', fontsize=20)
plt.plot(results[:,2], results[:,4]/(0.5*Uinf**2*Radius), 'r-', label=r'$C_t$')
plt.plot(results[:,2], results[:,3]/(0.5*Uinf**2*Radius), 'b-', label=r'$C_n$')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$C_t$ and $C_n$ [N]', fontsize=20)
plt.legend(fontsize=20)
#plt.savefig('Ct_Cn_lamb10.pdf', format='pdf', dpi=1000)

# Cq vs. r/R
fig1 = plt.figure(figsize=(12, 6))
plt.title(r'Torque, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R^2$', fontsize=20)
plt.plot(results[:,2], results[:,8]/(0.5*Uinf**2*Radius**2), 'g-', label=r'$C_q$')
plt.grid()
plt.xlabel('r/R', fontsize=20)
plt.ylabel(r'$C_q$ [N]', fontsize=20)
plt.legend(fontsize=20)
#plt.savefig('Cq_lamb10.pdf', format='pdf', dpi=1000)
