# import necessary libraries
from platypus import NSGAII, Problem, Real
import matplotlib.pyplot as plt
import numpy as np
import BEMmodel_optPR as BEMPR
import BEMmodel_optWT as BEMWT

""" ---- Code variables ---- """
N = 100 # annulus spacing

""" ---- Optimization variables ---- """
# Objective: For a given Nb and R, change the chord and the twist distribution 
# to maximise CP
WTOpt = False # Select what you want to optimize
PROpt = True

if WTOpt == True:
    
    CT = 0.75
    NB = 3
    Radius = 50
    Uinf = 10
    TSR = 8
    Omega = Uinf*TSR/Radius
    
    Variables = [CT, NB, Radius, Uinf, TSR, Omega]
    c_coefs = [0.63896187, -1.39013372,  0.5868371 ,  0.72921565, -0.86310821, 0.37711136]
    t_coefs = [-129.4313803 , 398.66429627, -484.78858014, 310.92637666, -125.77094024, 33.4854036]
    
    [CP, results, fnorm] = BEMWT.ExecuteBEM_optWT(PROpt, WTOpt, N, Variables, c_coefs, t_coefs)
    
elif PROpt == True:
    
    # Theoretical optimum
#    c_coefs = [0.63896187, -1.39013372,  0.5868371 ,  0.72921565, -0.86310821, 0.37711136]
#    t_coefs = [-129.4313803 , 398.66429627, -484.78858014, 310.92637666, -125.77094024, 33.4854036]
    def PropOptimized(vars):
        # Input data
        NB = 6
        Radius = 0.7
        Uinf = 60
        Omega = 1200/60*2*np.pi
        Params = [NB, Radius, Uinf, Omega]
        
        # Coefficients to optimize
        c1 = vars[0]
        c2 = vars[1] 
        c3 = vars[2]
        c4 = vars[3]
        c5 = vars[4]
        c6 = vars[5]
        t1 = vars[6]
        t2 = vars[7]
        t3 = vars[8]
        t4 = vars[9]
        t5 = vars[10]
        t6 = vars[11]
        
        c_coefs = [c1, c2, c3, c4, c5, c6]
        t_coefs = [t1, t2, t3, t4, t5, t6]
        
        [CP] = BEMPR.ExecuteBEM_optPR(N, Params, c_coefs, t_coefs)
        
        return[CP]

    # Creation of the optimisation algorithm
            
    problem = Problem(12, 1)
    problem.directions[:] = Problem.MAXIMIZE
    problem.types[:] = [Real(-2, 2), Real(-2, 2), Real(-2, 2), Real(-2, 2), Real(-2, 2), Real(-2, 2), 
                        Real(-1e3, 1e3), Real(-1e3, 1e3), Real(-1e3, 1e3), Real(-1e3, 1e3), Real(-1e3, 1e3), Real(-1e3, 1e3)]
    problem.function = PropOptimized
    
    algorithm = NSGAII(problem)
    algorithm.run(10^5)
    
    feasible_solutions = [s for s in algorithm.result if s.feasible]
    
    CP_it = [s.objectives[0] for s in algorithm.result]
    
    c1_it = [s.variables[0] for s in algorithm.result]
    c2_it = [s.variables[1] for s in algorithm.result]
    c3_it = [s.variables[2] for s in algorithm.result]
    c4_it = [s.variables[3] for s in algorithm.result]
    c5_it = [s.variables[4] for s in algorithm.result]
    c6_it = [s.variables[5] for s in algorithm.result]
    t1_it = [s.variables[6] for s in algorithm.result]
    t2_it = [s.variables[7] for s in algorithm.result]
    t3_it = [s.variables[8] for s in algorithm.result]
    t4_it = [s.variables[9] for s in algorithm.result]
    t5_it = [s.variables[10] for s in algorithm.result]
    t6_it = [s.variables[11] for s in algorithm.result]
    
    # Results after the optimisation
    m_CP = 0
    for i in range(len(CP_it)):
        # Also check that CP is lower than Betz limit !
        if (CP_it[i] > m_CP and CP_it[i] < 16/27):
            m_CP = CP_it[i]
            
            # Save the corresponding position of optimal twist and chord law
            pos = i
          
    display(m_CP)
    # The optimal chord and twist laws are
    c1_opt = c1_it[pos]
    c2_opt = c2_it[pos]
    c3_opt = c3_it[pos]
    c4_opt = c4_it[pos]
    c5_opt = c5_it[pos]
    c6_opt = c6_it[pos]
    t1_opt = t1_it[pos]
    t2_opt = t2_it[pos]
    t3_opt = t3_it[pos]
    t4_opt = t4_it[pos]
    t5_opt = t5_it[pos]
    t6_opt = t6_it[pos]
    
    delta_r_R = 1/N
    r_R = np.arange(0.25, 1+delta_r_R/2, delta_r_R)
    
    c_R_opt = c1_opt*r_R**5 + c2_opt*r_R**4 + c3_opt*r_R**3 +  \
                             c4_opt*r_R**2 + c5_opt*r_R + c6_opt
    twi_opt = t1_opt*r_R**5 + t2_opt*r_R**4 + t3_opt*r_R**3 +  \
                             t4_opt*r_R**2 + t5_opt*r_R + t6_opt
                             
    fig1 = plt.figure(figsize=(12, 6))
    plt.plot(r_R, c_R_opt, 'g--', label=r'c_R Optimised')
#    plt.plot(r_R, twi_opt, 'g--', label=r'$\theta$ Optimised')
    plt.xlabel(r'r_R', fontsize=20)
    plt.ylabel(r'c_R and $\theta$', fontsize=20)
    plt.grid()
    plt.legend(fontsize=20)

