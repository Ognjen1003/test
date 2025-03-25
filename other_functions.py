# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 04:37:30 2024

@author: domagoj
"""
import numpy as np

def Prandtl(Cp, mu, lambda_f):
    return Cp*mu/lambda_f

def dT(T0, Ts, Dp, Dt, G, Cp, L, lambda_i, lambda_f, lambda_s, Re_f, Prt, Prf, h):
    """
    Sukhov model for temperature drop in pipeline.
    
    Parameters
    ----------
    T0 : float
        temperature at distance 0 (inlet), K
    Ts : float
        temperature of soil, K
    Dp : float
        outer diameter of pipe, m
    Dt : float
        inner diameter of pipe, m   
    G : float
        mass flow rate, kg/s
    Cp : float
        heat capacity, kJ/kg·K
    L : float
        length of pipeline, km
    lambda_i : float
        thermal conductivity of isolation, W/m·K
    lambda_f : float
        thermal conductivity of fluid, W/m·K
    lambda_s : float
        thermal conductivity of soil, W/m·K
    Re_f : float
        Reynolds number
    Prt : float
        Prandtl number at tube temperature.
    Prf : float
        Prandtl number at fluid temperature.
    h : float
        depth of pipe underground, m

    Returns
    -------
    TL : float
        temperature at distance L (outlet).
    """
    alpha1 = 0.021*(lambda_f/Dt)*(Re_f**0.8)*(Prf**0.44)*(Prf/Prt)**0.25
    alpha2 = 2*lambda_s/((0.5*(Dt+Dp)*np.ln((2*h/(0.5*(Dt+Dp))+((2*h/(0.5*(Dt+Dp)))**2)-1))))
    K = (1/(1/(alpha1*Dt)+1/(alpha2*Dp)+np.ln(Dp/Dt)/(2*lambda_i)))*2/(Dt+Dp)
    TL = T0+(Ts-T0)*np.exp(-K*np.pi*Dp*L/(G*Cp))
    return TL

def d_Economic(Q, rho, mu):
    """
    Calculates optimal economic pipe diameter (Zhang, Z.X., Wang, G.X., 
            Massarotto, P., Rudolph, V., 2006. Optimization of pipeline
            transport for CO2 sequestration. 
            Energy Conversion and Management, 47: 702-715.)

    Parameters
    ----------
    Q : float
        mass flow rte kg/s.
    rho : float
        max density of the fluid in a pipe, kg/m3.
    mu : float
        max viscosity of the fluid in a pipe, Pas.
    Returns
    -------
    economic optimum diameter for CO2 transport, float
    """
    return (0.363*((Q/rho)**0.45)*(rho**0.13)*mu**0.025)

def wall_thickness(P_max, D, S = 483, F = 0.72 , E = 1):
    """
    Parameters
    ----------
    P_max : float
        Max. design pressure (MPa).
    D : float
        Internal diameter (m).
    S : float
        Pipe yield strength (MPa).
    F : float
        Safety (design) factor.
    E : float
        longitudilnal seam joint factor.
    Returns
    -------
    recommended wall thickness, mm
    """
    wt = 0.5*P_max*D/(S*F*E-P_max)
    return wt