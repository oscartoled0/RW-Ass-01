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
        a1=1-np.sqrt(CT1)/2;
        CT[a>a1] = CT1-4*(np.sqrt(CT1)-1)*(1-a[a>a1])
    
    return CT
     
def ainduction(CT):
    """
    This function calculates the induction factor 'a' as a function of thrust coefficient CT 
    including Glauert's correction
    """
    a = np.zeros(np.shape(CT))
    CT1=1.816;
    CT2=2*np.sqrt(CT1)-CT1
    a[CT>=CT2] = 1 + (CT[CT>=CT2]-CT1)/(4*(np.sqrt(CT1)-1))
    a[CT<CT2] = 0.5-0.5*np.sqrt(1-CT[CT<CT2])
    return a

def PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction):
    """
    This function calcualte steh combined tip and root Prandtl correction at agiven radial position 'r_R' (non-dimensioned by rotor radius), 
    given a root and tip radius (also non-dimensioned), a tip speed ratio TSR, the number lf blades NBlades and the axial induction factor
    """
    temp1 = -NBlades/2*(tipradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1-axial_induction)**2))
    Ftip = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Ftip[np.isnan(Ftip)] = 0
    temp1 = NBlades/2*(rootradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1-axial_induction)**2))
    Froot = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Froot[np.isnan(Froot)] = 0
    return Froot*Ftip, Ftip, Froot

# define function to determine load in the blade element
def loadBladeElement(vnorm, vtan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd):
    """
    calculates the load in the blade element
    """
    vmag2 = vnorm**2 + vtan**2
    inflowangle = np.arctan2(vnorm,vtan)
    alpha = inflowangle*180/np.pi + twist
    cl = np.interp(alpha, polar_alpha, polar_cl)
    cd = np.interp(alpha, polar_alpha, polar_cd)
    lift = 0.5*vmag2*cl*chord
    drag = 0.5*vmag2*cd*chord
    fnorm = lift*np.cos(inflowangle)+drag*np.sin(inflowangle)
    ftan = lift*np.sin(inflowangle)-drag*np.cos(inflowangle)
    gamma = 0.5*np.sqrt(vmag2)*cl*chord
    return fnorm , ftan, gamma, inflowangle*180/np.pi, alpha

#def solveStreamtube(Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, r_root, r_root2):
#    """
#    solve balance of momentum between blade element load and loading in the streamtube
#    input variables:
#    Uinf - wind speed at infinity
#    r1_R,r2_R - edges of blade element, in fraction of Radius ;
#    rootradius_R, tipradius_R - location of blade root and tip, in fraction of Radius ;
#    Radius is the rotor radius
#    Omega -rotational velocity
#    NBlades - number of blades in rotor
#    """
#    Area = np.pi*((r2_R*Radius)**2-(r1_R*Radius)**2) #  area streamtube
#    r_R = (r1_R+r2_R)/2 # centroide
#    # initiatlize variables
#    a = 0.3 # axial induction
#    aline = 0.0 # tangential induction factor
#    
#    Niterations = 100
#    count = 1
#    Erroriterations = 1e-5# error limit for iteration process, in absolute value of induction
#    conv = False
#    while conv== False and count < Niterations:
#        # ///////////////////////////////////////////////////////////////////////
#        # // this is the block "Calculate velocity and loads at blade element"
#        # ///////////////////////////////////////////////////////////////////////
#        Urotor = Uinf*(1-a) # axial velocity at rotor
#        Utan = (1+aline)*Omega*r_R*Radius # tangential velocity at rotor
#        # calculate loads in blade segment in 2D (N/m)
#        fnorm, ftan, gamma, inflowangle, alpha = loadBladeElement(Urotor, Utan, r_R,chord, twist, polar_alpha, polar_cl, polar_cd)
#        load3Daxial =fnorm*Radius*(r2_R-r1_R)*NBlades # 3D force in axial direction
#        # load3Dtan =loads[1]*Radius*(r2_R-r1_R)*NBlades # 3D force in azimuthal/tangential direction (not used here)
#      
#        # ///////////////////////////////////////////////////////////////////////
#        # //the block "Calculate velocity and loads at blade element" is done
#        # ///////////////////////////////////////////////////////////////////////
#
#        # ///////////////////////////////////////////////////////////////////////
#        # // this is the block "Calculate new estimate of axial and azimuthal induction"
#        # ///////////////////////////////////////////////////////////////////////
#        # // calculate thrust coefficient at the streamtube 
#        CT = load3Daxial/(0.5*Area*Uinf**2)
#        
#        # calculate new axial induction, accounting for Glauert's correction
#        anew =  ainduction(CT)
#        
#        # correct new axial induction with Prandtl's correction
#        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, anew);
#        if (Prandtl < 1e-4): 
#            Prandtl = 1e-4 # avoid divide by zero
#        anew = anew/Prandtl # correct estimate of axial induction
#        a = 0.75*a+0.25*anew # for improving convergence, weigh current and previous iteration of axial induction
#          
#         # calculate aximuthal induction
#        aline = ftan*NBlades/(2*np.pi*Uinf*(1-a)*Omega*2*(r_R*Radius)**2)
#        alin_new =aline/Prandtl # correct estimate of azimuthal induction with Prandtl's correction
#        
#        aline = 0.75*alin_new + 0.25*aline
#        # ///////////////////////////////////////////////////////////////////////////
#        # // end of the block "Calculate new estimate of axial and azimuthal induction"
#        # ///////////////////////////////////////////////////////////////////////
#        
#        # Avoid stack overflows
#        if a > 0.95:
#            a = 0.95
#        elif a < 0:
#            a = 0
#            
#        if aline > 0.95:
#            aline = 0.95
#        elif aline < 0:
#            aline = 0
#            
#        # if r_R == (r_root+r_root2)/2:
#        #     display('a/at/Pr')
#        #     display(a)
#        #     display(aline)
#        #     display(fnorm)
#            
#        count = count + 1
#        #// test convergence of solution, by checking convergence of axial induction
#        if (np.abs(a-anew) < Erroriterations): 
#            conv = True
#             
#        
#        fQ = ftan*r_R
#    return [a , aline, r_R, fnorm , ftan, gamma, inflowangle, alpha, fQ]

def solveStreamtubeWT(Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, r_root, r_root2, CT_r):
    """
    Solves to keep C_T = 0.75
    """
    Area = np.pi*((r2_R*Radius)**2-(r1_R*Radius)**2) #  area streamtube
    r_R = (r1_R+r2_R)/2 # centroide
    a =  ainduction(CT_r)
    fnorm = CT_r*Area*0.5*Uinf**2
    vnorm = Uinf*(1-a) # axial velocity at rotor
    aline = a*(1-a)/(TSR**2*r_R**2)

#    Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, a);
#    if (Prandtl < 1e-4): 
#        Prandtl = 1e-4 # avoid divide by zero
#    a = a/Prandtl # correct estimate of axial induction
    
    aline = 0.0 # Guess aline
    Niterations = 10
    count = 1
    conv = False
    while conv== False and count<Niterations:
        vtan = (1+aline)*Omega*r_R*Radius # tangential velocity at rotor
        vmag2 = vnorm**2 + vtan**2
        inflowangle = np.arctan2(vnorm,vtan)
        alpha = inflowangle*180/np.pi + twist
        cl = np.interp(alpha, polar_alpha, polar_cl)
        cd = np.interp(alpha, polar_alpha, polar_cd)
        lift = 0.5*vmag2*cl*chord
        drag = 0.5*vmag2*cd*chord
        fnorm2 = lift*np.cos(inflowangle)+drag*np.sin(inflowangle)
        ftan = lift*np.sin(inflowangle)-drag*np.cos(inflowangle)
        gamma = 0.5*np.sqrt(vmag2)*cl*chord
        fQ = ftan*r_R
        
        if abs(fnorm2 - fnorm)<0.001:
            conv == True
        else:
            aline = ftan*NBlades/(2*np.pi*Uinf*(1-a)*Omega*2*(r_R*Radius)**2)
                        
        # if aline > 0.95:
        #     aline = 0.95
        # elif aline < 0:
        #     aline = 0
            
        count = count + 1
        display('abs(fnorm2 - fnorm)')
        display(abs(fnorm2 - fnorm))

    return [a , aline, r_R, fnorm2 , ftan, gamma, inflowangle, alpha, fQ]

""" ------- Sections ------- """
def ExecuteBEM_optWT(PROpt, WTOpt, N, Variables, c, t):

    """ POLARS """
        
    if (WTOpt == True and PROpt == False):
        # Wind turbine rotor
        airfoil = 'DU95W180.cvs'
        data1=pd.read_csv(airfoil, header=0,
                            names = ["alfa", "cl", "cd", "cm"],  sep='\s+')
        polar_alpha = data1['alfa'][:]
        polar_cl = data1['cl'][:]
        polar_cd = data1['cd'][:]
    elif (PROpt == True and WTOpt == False):
        # Propeller case
        airfoil = 'ARAD8polar.csv'
        data1=pd.read_csv(airfoil, header=0,
                          names = ["alfa", "cl", "cd", "cm"],  sep='\s+')
        polar_alpha = np.flipud(-data1['alfa'][:])
        polar_cl = np.flipud(-data1['cl'][:])
        polar_cd = np.flipud(data1['cd'][:])

        
    """ BLADE GEOMETRY AND FLOW CONDITIONS"""
    
    delta_r_R = 1/N
    if PROpt == True:
        
        NBlades = Variables[0], Radius = Variables[1], Uinf = Variables[2], Omega = Variables[3]
        r_R = np.arange(0.25, 1+delta_r_R/2, delta_r_R)
        chord_distribution = c[0]*r_R**5 + c[1]*r_R**4 + c[2]*r_R**3 + c[3]*r_R**3 + c[2]*r_R**2 + c[1]*r_R + c[0]
        twist_distribution = t[0]*r_R**5 + t[1]*r_R**4 + t[2]*r_R**3 + t[3]*r_R**3 + t[2]*r_R**2 + t[1]*r_R + t[0]
        
    elif WTOpt == True:
        
        CT = Variables[0]
        NBlades = Variables[1]
        Radius = Variables[2]
        Uinf = Variables[3]
        TSR = Variables[4]
        Omega = Variables[5]
        r_R = np.arange(0.2, 1+delta_r_R/2, delta_r_R)
        chord_distribution = c[0]*r_R**5 + c[1]*r_R**4 + c[2]*r_R**3 + c[3]*r_R**3 + c[2]*r_R**2 + c[1]*r_R + c[0]
        twist_distribution = t[0]*r_R**5 + t[1]*r_R**4 + t[2]*r_R**3 + t[3]*r_R**3 + t[2]*r_R**2 + t[1]*r_R + t[0]
    
    """ SOLVE BEM """
            
    areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*Radius**2
    dr = (r_R[1:]-r_R[:-1])*Radius
    
    CT_r = np.zeros(len(r_R)-1)
    for i in range(len(r_R)-1):
        r_R_med = (r_R[i]+r_R[i+1])/2
        m = CT/(1/2 - 0.2)
        CT_r[i] = m*(r_R_med - 0.2) # That is just a linear curve for the CT
        chord = np.interp((r_R[i]+r_R[i+1])/2, r_R, chord_distribution)
        twist = np.interp((r_R[i]+r_R[i+1])/2, r_R, twist_distribution)
        
        display('New r/R!')
        results =np.zeros([len(r_R)-1,9]) 
        results[i,:] = solveStreamtubeWT(Uinf, r_R[i], r_R[i+1], r_R[0], r_R[-1] , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, r_R[0], r_R[1], CT_r[i])
        
        CP = np.sum(dr*results[:,4]*results[:,2]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))   
    
    return[CP, results, fnorm]
        