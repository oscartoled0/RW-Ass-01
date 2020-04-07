# import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as m

""" ------- Functions ------- """

def CTfunction(a, glauert = False):
    """
    This function calculates the thrust coefficient as a function of induction factor 'a'
    'glauert' defines if the Glauert correction for heavily loaded rotors should be used; default value is false
    """
    CT = np.zeros(np.shape(a))
    CT = 4*a*(1-a)  
    if glauert:
        CT1=1.816;
        a1=1-np.sqrt(CT1)/2
        if a>a1:
            CT = CT1-4*(np.sqrt(CT1)-1)*(1-a)
    
    return CT

def ainduction(CT):
    """
    This function calculates the induction factor 'a' as a function of thrust coefficient CT 
    including Glauert's correction
    """
    a = np.zeros(np.shape(CT))
    CT1=1.816;
    CT2=2*np.sqrt(CT1)-CT1
    if CT>=CT2:
        a = - 1 + (CT1-CT)/(4*(np.sqrt(CT1)-1))
    elif CT<CT2:
        a = 0.5 - 0.5*np.sqrt(1-CT)
        
    return a

def PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction):
    """
    This function calcualte steh combined tip and root Prandtl correction at agiven radial position 'r_R' (non-dimensioned by rotor radius), 
    given a root and tip radius (also non-dimensioned), a tip speed ratio TSR, the number lf blades NBlades and the axial induction factor
    """
    temp1 = -NBlades/2*(tipradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1+axial_induction)**2))
    Ftip = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Ftip[np.isnan(Ftip)] = 0
    temp1 = NBlades/2*(rootradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1+axial_induction)**2))
    Froot = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Froot[np.isnan(Froot)] = 0
    return Froot*Ftip, Ftip, Froot

""" ------- Sections ------- """
def ExecuteBEM_opt(N, Variables):

    """ POLARS """
        
    # Propeller case
    airfoil = 'ARAD8polar.csv'
    data1=pd.read_csv(airfoil, header=0,
                      names = ["alfa", "cl", "cd", "cm"],  sep='\s+')
    polar_alpha = data1['alfa'][:]
    polar_cl = data1['cl'][:]
    polar_cd = data1['cd'][:]
        
    """ MAXIMUM EFFICIENCY """
    
    polar_E = np.zeros(polar_alpha.size)
    max_E = 0
    for i in range(polar_alpha.size):
        polar_E = polar_cl[i]/polar_cd[i]
        if (max_E < polar_E) and abs(polar_alpha[i])<30:
            max_E = polar_E
            alpha_E = polar_alpha[i]
            cl_E = polar_cl[i]
            cd_E = polar_cd[i]
            
    # display(alpha_E)

    """ BLADE GEOMETRY AND FLOW CONDITIONS"""
    
    delta_r_R = 1/N
    
    CT = Variables[0]
    NBlades = Variables[1]
    Radius = Variables[2]
    Uinf = Variables[3]
    TSR = Variables[4]
    Omega = Variables[5]
    
    r = np.arange(0.25, 1+delta_r_R/2, delta_r_R)
        
    """ CONSTRUCT OPTIMIZED BLADE """
    
    # We define the vectors of parameters
    a = np.zeros(r.size-1)
    aline = np.zeros(r.size-1)
    phi = np.zeros(r.size-1)
    twist = np.zeros(r.size-1)
    chord = np.zeros(r.size-1)
    fnorm = np.zeros(r.size-1)
    ftan = np.zeros(r.size-1)
    gamma = np.zeros(r.size-1)
    r_R = np.zeros(r.size-1)
    vmag2 = np.zeros(r.size-1)
    
    for i in range(r.size-1):
        
        r_R[i] = (r[i]+r[i+1])/2
        Area = np.pi*((r[i+1]*Radius)**2-(r[i]*Radius)**2)
        
        a[i] = ainduction(CT)
        # display(a[i])
            
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(r_R[i], r[0], r[-1], Omega*Radius/Uinf, NBlades, a[i]);
        if (Prandtl < 1e-4): 
            Prandtl = 1e-4 # avoid divide by zero
        a[i] = a[i]/Prandtl
        
        if a[i]>0.95: # a cannot be larger than 0.95 because if not reverse flow
            a[i] = 0.95
        elif a[i] < 0:
            a[i] = 0
            
        CTd=CTfunction(a[i], glauert = False)
        
        aline[i] = a[i]*(1+a[i])/(TSR**2*r_R[i]**2)
        vnorm = Uinf*(1+a[i])
        vtan = (1-aline[i])*Omega*r_R[i]*Radius # tangential velocity at rotor
        phi[i] = np.arctan2(vnorm,vtan)
        
        # Now the optimum twist can be derived from the inflow angle for a_opt and alpha for E_max:
        twist[i] = alpha_E + phi[i]*180/np.pi
        
        # Now we calculate the aerodynamic forces in the blade in order to derive the optimum chord
        fnorm[i] = CTd*(0.5*Area*Uinf**2)/(Radius*(r[i+1]-r[i])*NBlades)
        vmag2[i] = vnorm**2 + vtan**2
        chord[i] = fnorm[i]/(0.5*vmag2[i]*(cl_E*m.cos(phi[i]) + cd_E*m.sin(phi[i])))
        ftan[i] = 0.5*vmag2[i]*chord[i]*(cl_E*np.sin(phi[i])-cd_E*np.cos(phi[i]))
        gamma[i] = 0.5*np.sqrt(vmag2[i])*cl_E*chord[i]
        # gamma[i] = (4*m.pi*a[i]*(1-a[i])*Uinf**2)/Omega
        
    dr = (r[1:]-r[:-1])*Radius
    CP = np.sum(dr*ftan*r_R*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))   
    CT = np.sum(dr*fnorm*NBlades/(0.5*Uinf**2*np.pi*Radius**2))

    return[CP, CT, fnorm, ftan, gamma, chord, twist, r_R, phi, a, aline, vmag2, r]