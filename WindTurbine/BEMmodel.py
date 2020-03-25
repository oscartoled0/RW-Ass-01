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

def solveStreamtube(Uinf, r1_R, r2_R, rootradius_R, tipradius_R , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, r_root, r_root2):
    """
    solve balance of momentum between blade element load and loading in the streamtube
    input variables:
    Uinf - wind speed at infinity
    r1_R,r2_R - edges of blade element, in fraction of Radius ;
    rootradius_R, tipradius_R - location of blade root and tip, in fraction of Radius ;
    Radius is the rotor radius
    Omega -rotational velocity
    NBlades - number of blades in rotor
    """
    Area = np.pi*((r2_R*Radius)**2-(r1_R*Radius)**2) #  area streamtube
    r_R = (r1_R+r2_R)/2 # centroide
    # initiatlize variables
    a = 0.6 # axial induction
    aline = 0.0 # tangential induction factor
    
    Niterations = 10^5
    count = 1
    Erroriterations = 1e-5# error limit for iteration process, in absolute value of induction
    conv = False
    while conv== False and count < Niterations:
        # ///////////////////////////////////////////////////////////////////////
        # // this is the block "Calculate velocity and loads at blade element"
        # ///////////////////////////////////////////////////////////////////////
        Urotor = Uinf*(1-a) # axial velocity at rotor
        Utan = (1+aline)*Omega*r_R*Radius # tangential velocity at rotor
        # calculate loads in blade segment in 2D (N/m)
        fnorm, ftan, gamma, inflowangle, alpha = loadBladeElement(Urotor, Utan, r_R,chord, twist, polar_alpha, polar_cl, polar_cd)
        load3Daxial =fnorm*Radius*(r2_R-r1_R)*NBlades # 3D force in axial direction
        # load3Dtan =loads[1]*Radius*(r2_R-r1_R)*NBlades # 3D force in azimuthal/tangential direction (not used here)
      
        # ///////////////////////////////////////////////////////////////////////
        # //the block "Calculate velocity and loads at blade element" is done
        # ///////////////////////////////////////////////////////////////////////

        # ///////////////////////////////////////////////////////////////////////
        # // this is the block "Calculate new estimate of axial and azimuthal induction"
        # ///////////////////////////////////////////////////////////////////////
        # // calculate thrust coefficient at the streamtube 
        CT = load3Daxial/(0.5*Area*Uinf**2)
        
        # calculate new axial induction, accounting for Glauert's correction
        anew =  ainduction(CT)
        
        # correct new axial induction with Prandtl's correction
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, Omega*Radius/Uinf, NBlades, anew);
        if (Prandtl < 1e-4): 
            Prandtl = 1e-4 # avoid divide by zero
        anew = anew/Prandtl # correct estimate of axial induction
        a = 0.75*a+0.25*anew # for improving convergence, weigh current and previous iteration of axial induction
          
         # calculate aximuthal induction
        aline = ftan*NBlades/(2*np.pi*Uinf*(1-a)*Omega*2*(r_R*Radius)**2)
        alin_new =aline/Prandtl # correct estimate of azimuthal induction with Prandtl's correction
        
        aline = 0.75*alin_new + 0.25*aline
        # ///////////////////////////////////////////////////////////////////////////
        # // end of the block "Calculate new estimate of axial and azimuthal induction"
        # ///////////////////////////////////////////////////////////////////////
        
        # Avoid stack overflows
        if a > 0.95:
            a = 0.95
        elif a < 0:
            a = 0
            
        if aline > 0.95:
            aline = 0.95
        elif aline < 0:
            aline = 0
            
        # if r_R == (r_root+r_root2)/2:
        #     display('a/at/Pr')
        #     display(a)
        #     display(aline)
        #     display(fnorm)
            
        count = count + 1
        #// test convergence of solution, by checking convergence of axial induction
        if (np.abs(a-anew) < Erroriterations): 
            conv = True             
        
        fQ = ftan*r_R
    return [a , aline, r_R, fnorm , ftan, gamma, inflowangle, alpha, fQ]

""" ------- Sections ------- """
def ExecuteBEM(PR, WT, N, plotter, mode = 'constant'):

        """ 1. plot CT as a function of induction "a", with and without Glauert correction """
        # define a as a range
        a = np.arange(-.5,1,.01)
        CTmom = CTfunction(a) # CT without correction
        CTglauert = CTfunction(a, True) # CT with Glauert's correction
        a2 = ainduction(CTglauert)
        
        if plotter == True:
            fig1 = plt.figure(figsize=(12, 6))
            plt.plot(a, CTmom, 'k-', label='$C_T$')
            plt.plot(a, CTglauert, 'b--', label='$C_T$ Glauert')
            plt.plot(a, CTglauert*(1-a), 'g--', label='$C_P$ Glauert')
            plt.xlabel('a', fontsize=20)
            plt.ylabel(r'$C_T$ and $C_P$', fontsize=20)
            plt.grid()
            plt.legend(fontsize=20)
            plt.savefig('Glauert.pdf', format='pdf', dpi=1000)
        
        """ 2. plot Prandtl tip, root and combined correction for a number of blades and induction 'a', over the non-dimensioned radius """
        r_R = np.arange(0.1, 1, .01)
        a = np.zeros(np.shape(r_R))+0.3
        Prandtl, Prandtltip, Prandtlroot = PrandtlTipRootCorrection(r_R, 0.1, 1, 7, 3, a)
        
        if plotter == True:
            fig1 = plt.figure(figsize=(12, 6))
            plt.plot(r_R, Prandtl, 'r-', label='Prandtl')
            plt.plot(r_R, Prandtltip, 'g.', label='Prandtl tip')
            plt.plot(r_R, Prandtlroot, 'b.', label='Prandtl root')
            plt.xlabel('r/R', fontsize=20)
            plt.legend(fontsize=20)
            plt.grid()
            plt.savefig('Prandtl.pdf', format='pdf', dpi=1000)

        """ 3. import polar """
        
        if (WT == True and PR == False):
            # Wind turbine rotor
            airfoil = 'DU95W180.cvs'
            data1=pd.read_csv(airfoil, header=0,
                                names = ["alfa", "cl", "cd", "cm"],  sep='\s+')
            polar_alpha = data1['alfa'][:]
            polar_cl = data1['cl'][:]
            polar_cd = data1['cd'][:]
        elif (PR == True and WT == False):
            # Propeller case
            airfoil = 'ARAD8polar.csv'
            data1=pd.read_csv(airfoil, header=0,
                              names = ["alfa", "cl", "cd", "cm"],  sep='\s+')
#            polar_alpha = data1['alfa'][:]
#            polar_cl = data1['cl'][:]
#            polar_cd = data1['cd'][:]
            polar_alpha = np.flipud(-data1['alfa'][:])
            polar_cl = np.flipud(-data1['cl'][:])
            polar_cd = np.flipud(data1['cd'][:])

        # plot polars of the airfoil C-alfa and Cl-Cd
        if plotter == True:
            fig, axs = plt.subplots(1, 2, figsize=(12, 6))
            axs[0].plot(polar_alpha, polar_cl)
            axs[0].set_xlim([-30,30])
            axs[0].set_xlabel(r'$\alpha$',fontsize=20)
            axs[0].set_ylabel(r'$C_l$',fontsize=20)
            axs[0].grid()
            axs[1].plot(polar_cd, polar_cl)
            axs[1].set_xlim([0,.1])
            axs[1].set_xlabel(r'$C_d$',fontsize=20)
            axs[1].grid()
            plt.savefig('Polar_prop.pdf', format='pdf', dpi=1000)

            
        """ 4. Define the blade geometry """
        if (WT == True and PR == False):
            # Wind turbine rotor
            # Basic rotor specs
            Radius = 50
            NBlades = 3
            
            # Blade specs
            if mode == 'constant':
                delta_r_R = 1/N
                r_R = np.arange(0.2, 1+delta_r_R/2, delta_r_R)
            elif mode == 'cosinus': 
                r_R = np.zeros(N)
                for i in range(N):
                    r_R[i] = 0.2 + (0.4)*(1-m.cos(((i)/(N-1))*m.pi))
                    
                # plt.figure(figsize=(12, 6))
                # plt.scatter(r_R,r_R)
                # plt.title('Radial distribution')
            else:
                display('Please enter a valid distribution mode.')

            pitch = 2 # degrees
            twist_distribution = -14*(1-r_R)+pitch # degrees
            chord_distribution = 3*(1-r_R)+1 # meters
            
            # Operational specs
            Uinf = 10 # unperturbed wind speed in m/s
            TSR = 10 # tip speed ratio
            Omega = Uinf*TSR/Radius
        elif (PR == True and WT == False):
            # Propeller rotor
            # Basic rotor specs
            Radius = 0.7
            NBlades = 6
            
            # Blade specs
            if mode == 'constant':
                delta_r_R = 1/N
                r_R = np.arange(0.25, 1+delta_r_R/2, delta_r_R)
            elif mode == 'cosinus': 
                r_R = np.zeros(N)
                for i in range(N):
                    r_R[i] = 0.25 + (0.75/2)*(1-m.cos(((i)/(N-1))*m.pi))
                    
                # plt.figure(figsize=(12, 6))
                # plt.scatter(r_R,r_R)
                # plt.title('Radial distribution')
            else:
                display('Please enter a valid distribution mode.')
                
            twist_distribution = (-50*(r_R)+35+46) # degrees
            chord_distribution = 0.18-0.06*(r_R) # meters
            
            # Operational specs
            Uinf = 60 # unperturbed wind speed in m/s
            Omega = 1200/60*2*np.pi
            TSR = Omega*Radius/Uinf
            # h = 2000 #m
        # solve BEM model
        
        results =np.zeros([len(r_R)-1,9]) 
        
        for i in range(len(r_R)-1):
            chord = np.interp((r_R[i]+r_R[i+1])/2, r_R, chord_distribution)
            twist = np.interp((r_R[i]+r_R[i+1])/2, r_R, twist_distribution)
            
            results[i,:] = solveStreamtube(Uinf, r_R[i], r_R[i+1], r_R[0], r_R[-1] , Omega, Radius, NBlades, chord, twist, polar_alpha, polar_cl, polar_cd, r_R[0], r_R[1] )
            
        """ 5. Plot results """
        
        # # To avoid INF in forces calculation
        # for i in range(results[:,3].size):
        #     if results[i,3] > 0.5:
        #         results[i,3] = 0
        #         results[i,4] = 0
        #         results[i,5] = 0
            
        areas = (r_R[1:]**2-r_R[:-1]**2)*np.pi*Radius**2
        dr = (r_R[1:]-r_R[:-1])*Radius
        CT = np.sum(dr*results[:,3]*NBlades/(0.5*Uinf**2*np.pi*Radius**2))
        CP = np.sum(dr*results[:,4]*results[:,2]*NBlades*Radius*Omega/(0.5*Uinf**3*np.pi*Radius**2))
            
        
        if plotter == True:
            print("CT is ", CT)
            print("CP is ", CP)    
            
            fig1 = plt.figure(figsize=(12, 6))
            plt.title('Axial and tangential induction')
            plt.plot(results[:,2], results[:,0], 'r-', label=r'$a$')
            plt.plot(results[:,2], results[:,1], 'g--', label=r'$a^,$')
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()
            
            fig1 = plt.figure(figsize=(12, 6))
            plt.title(r'Normal and tagential force, non-dimensioned by $\frac{1}{2} \rho U_\infty^2 R$')
            plt.plot(results[:,2], results[:,3]/(0.5*Uinf**2*Radius), 'r-', label=r'Fnorm')
            plt.plot(results[:,2], results[:,4]/(0.5*Uinf**2*Radius), 'g--', label=r'Ftan')
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()
            
            
            fig1 = plt.figure(figsize=(12, 6))
            plt.title(r'Circulation distribution, non-dimensioned by $\frac{\pi U_\infty^2}{\Omega * NBlades } $')
            plt.plot(results[:,2], results[:,5]/(np.pi*Uinf**2/(NBlades*Omega)), 'r-', label=r'$\Gamma$')
            plt.grid()
            plt.xlabel('r/R')
            plt.legend()
            plt.show()
        
        return[results, CT, CP, Uinf, Radius]
        