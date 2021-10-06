import numpy as np

def wave_power(H,T):
    '''Deepwater wave power'''
    rho = 1030. #kg/m3
    g = 9.81 # m/s/s
    n = 0.5 # deepwater
    L = g*T**2 / (2.*np.pi)
    E = rho*g*H**2*L/8.
    C = L/T
    P = E*C*n # Watts/m
    return P, L, C

def nondim_settling_vel(H,T,ws=0.05):
    '''Dean, 1973; Gourlay, 1968'''
    ws_nd = H/(ws*T)
    return(ws_nd)

def calcR2(H,T,slope,igflag=0):
    """
    %
    % [R2,S,setup, Sinc, SIG, ir] = calcR2(H,T,slope,igflag);
    %
    % Calculated 2% runup (R2), swash (S), setup (setup), incident swash (Sinc)
    % and infragravity swash (SIG) elevations based on parameterizations from runup paper
    % also Iribarren (ir)
    % August 2010 - Included 15% runup (R16) statistic that, for a Guassian distribution,
    % represents mean+sigma. It is calculated as R16 = setup + swash/4.
    % In a wave tank, Palmsten et al (2010) found this statistic represented initiation of dune erosion.
    %
    %
    % H = significant wave height, reverse shoaled to deep water
    % T = deep-water peak wave period
    % slope = radians
    % igflag = 0 (default)use full equation for all data
    %        = 1  use dissipative-specific calculations when dissipative conditions exist (Iribarren < 0.3)
    %        = 2  use dissipative-specific (IG energy) calculation for all data
    %
    % based on:
    %  Stockdon, H. F., R. A. Holman, P. A. Howd, and J. Sallenger A. H. (2006),
    %    Empirical parameterization of setup, swash, and runup,
    %    Coastal Engineering, 53, 573-588.
    % author: hstockdon@usgs.gov
    # Converted to Python by csherwood@usgs.gov
    """
    g = 9.81

    # make slopes positive!
    slope = np.abs(slope)

    # compute wavelength and Iribarren
    L = (g*T**2) / (2.*np.pi)
    sqHL = np.sqrt(H*L)
    ir = slope/np.sqrt(H/L)

    if igflag == 2:                     # use dissipative equations (IG) for ALL data
        R2 = 1.1*(0.039 * sqHL)
        S = 0.046*sqHL
        setup = 0.016*sqHL

    elif igflag == 1 and ir < 0.3:      # if dissipative site use diss equations
        R2 = 1.1*(0.039 * sqHL)
        S = 0.046*sqHL
        setup = 0.016*sqHL

    else:                               # if int/ref site, use full equations
        setup = 0.35*slope*sqHL
        Sinc = 0.75*slope*sqHL
        SIG = 0.06*sqHL
        S = np.sqrt(Sinc**2 + SIG**2)
        R2 = 1.1*(setup + S/2.)
        R16 = 1.1*(setup + S/4.)

    return R2, S, setup, Sinc, SIG, ir, R16
