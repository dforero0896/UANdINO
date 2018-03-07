#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#Constants
#Mass differences and mixing angles
kuoPantaleone=0
ohlsson=0
real = 1
if kuoPantaleone:
    dM32 = 1e-4
    dm21 = 1e-8
    theta1 = np.deg2rad(45.)
    theta2 = np.deg2rad(5.)
    theta3 = np.deg2rad(45.)
elif ohlsson:
    dM32 = 3.2e-3
    dm21 = 0.
    theta1 = np.deg2rad(45.)
    theta2 = np.deg2rad(5.)
    theta3 = np.deg2rad(45.)
elif real:
    dM32 = 2.45e-3
    dm21=7.53e-5
    theta1 = 0.78539816339744839
    theta2 = 0.1454258194533693
    theta3 = 0.5872523687443223;
# U matrix Elements
Ue1 = np.cos(theta2)*np.cos(theta3)
Ue2 = np.sin(theta3)*np.cos(theta2)
Ue3 = np.sin(theta2)
Umu1=-np.sin(theta3)*np.cos(theta1)-np.sin(theta1)*np.sin(theta2)*np.cos(theta3)
Umu2=np.cos(theta1)*np.cos(theta3)-np.sin(theta1)*np.sin(theta2)*np.sin(theta3)
Umu3=np.sin(theta1)*np.cos(theta2)
Ut1=np.sin(theta1)*np.sin(theta3)-np.sin(theta2)*np.cos(theta1)*np.cos(theta3)
Ut2=-np.sin(theta1)*np.cos(theta3)-np.sin(theta2)*np.sin(theta3)*np.cos(theta1)
Ut3=np.cos(theta1)*np.cos(theta2)
CKM = np.matrix([[Ue1, Ue2, Ue3], [Umu1, Umu2, Umu3], [Ut1, Ut2, Ut3]])

def sun_rho(r):
  return (200.)*np.exp(-np.abs(r)/66000); #g/cm^3


def sun_density(r):
    '''Returns density of the sun when given a coordinate r from it's center.'''
    x = 1.972e-16
    n0 = 245*6.022e23
    r0 = 6.57e5/10.54
    Gfoverhc3 = 1.1663787e-5
    ne = n0*np.exp(-r/r0)
    return np.sqrt(2.)*Gfoverhc3*ne*x*x*x*1e6*1e9

def fig_1_density(r):
  '''Returns the density of the Earth given a coordinate r from it's center following the convention provided in fig_1 of Ohlsonn's paper.'''
  dist = np.abs(r-(-6371))
  if(dist < 2885.):
      return 1.7e-13
  elif(dist>=2885. and dist<=2885.+6972.):
      return 4.35e-13
  else:
      return 1.7e-13

def fig_2_density(r):
  '''Returns the density of the Earth given a coordinate r from it's center following the convention provided in fig_2 of Ohlsonn's paper.'''
  dist = np.abs(r-(-6371));
  if(dist < 2885.):
      return 1.7e-13
  elif(dist>=2885. and dist<=2885.+6972.):
      return  2.1e-13
  else:
      return 1.7e-13

def fig_3_density(r):
  '''Returns the density of the Earth given a coordinate r from it's center following the convention provided in fig_3 of Ohlsonn's paper.'''
  dist = np.abs(r-(-6371));
  if(dist < 2885.):
      return 3.8e-14
  elif(dist>=2885. and dist<=2885.+6972.):
      return  7.6e-14
  else:
      return 3.8e-14

def fig_4_density(r):
  '''Returns the density of the Earth given a coordinate r from it's center following the convention provided in fig_4 of Ohlsonn's paper.'''
  return 3e-13

def fig_5_density(r):
  '''Returns the density of the Earth given a coordinate r from it's center following the convention provided in fig_5 of Ohlsonn's paper.'''
  return 1.7e-13

def fig_6_density(r):
  '''Returns the density of the Earth given a coordinate r from it's center following the convention provided in fig_6 of Ohlsonn's paper.'''
  dist = np.abs(r-(-6371));
  return 3.8e-13*(1e-3 +dist/12742.)

def density_to_potential(dty, antineutrino):
    to_return = (1./np.sqrt(2))*dty*1e-3*8.96189e-47*1e9   /1.672e-27
    if antineutrino:
        return -1*to_return
    else:
        return to_return

def longitude_units_conversion(lon_in_km):
    '''Transforms distances in km'''
    return lon_in_km*1e3/(1.972e-7)


def calculateOperator(neutrinoEnergy, A, L):
    '''Calculates the time-evolution operator for a(n) (anti)neutrino with energy "neutrinoEnergy" passing through a potential "A" and a distance "L".'''
    E21 = dm21/(2*neutrinoEnergy)
    E32 = dM32/(2*neutrinoEnergy)
    E12=-E21;
    E23=-E32;
    E31=E12-E23;
    E13=-E31;
    #Elements of the Tmatrix in mass basis
    T_11=A*Ue1*Ue1-(1./3)*A+(1./3)*(E12+E13);
    T_12=A*Ue1*Ue2;
    T_13=A*Ue1*Ue3;
    T_21=T_12;
    T_22=A*Ue2*Ue2-(1./3)*A+(1./3)*(E21+E23);
    T_23=A*Ue2*Ue3;
    T_31=T_13;
    T_32=T_23;
    T_33=A*Ue3*Ue3-(1./3)*A+(1./3)*(E31+E32);
    T_mass_mat = np.matrix([[T_11, T_12, T_13], [T_21, T_22, T_23], [T_31, T_32, T_33]])
    T_sq_11=(1./3)*(A*A*(Ue1*Ue1+(1./3))+2*A*(Ue1*Ue1-(1./3))*(E12+E13)+(1./3)*(E12+E13)*(E12+E13))
    T_sq_12=(1./3)*Ue1*Ue2*A*(A+E13+E23)
    T_sq_13=(1./3)*Ue1*Ue3*A*(A+E12+E32)
    T_sq_21=T_sq_12
    T_sq_22=(1./3)*(A*A*(Ue2*Ue2+(1./3))+2*A*(Ue2*Ue2-(1./3))*(E21+E23)+(1./3)*(E21+E23)*(E21+E23))
    T_sq_23=(1./3)*Ue2*Ue3*A*(A+E21+E31)
    T_sq_31=T_sq_13
    T_sq_32=T_sq_23
    T_sq_33=(1./3)*(A*A*(Ue3*Ue3+(1./3))+2*A*(Ue3*Ue3-(1./3))*(E31+E32)+(1./3)*(E31+E32)*(E31+E32))
    T_sq_mass_mat = np.matrix([[T_sq_11, T_sq_12, T_sq_13], [T_sq_21, T_sq_22, T_sq_23], [T_sq_31, T_sq_32, T_sq_33]])

    T_flav_mat = np.matmul(np.matmul(CKM,T_mass_mat),CKM.H)
    T_sq_flav_mat = np.matmul(np.matmul(CKM,T_sq_mass_mat),CKM.H)

    #Calculate c's
    c1 = -A*A/3 + (A/(6*neutrinoEnergy))*(Ue1*Ue1*(dM32+2*dm21)+Ue2*Ue2*(dM32-dm21)-Ue3*Ue3*(2*dM32+dm21))-(1./(12*neutrinoEnergy*neutrinoEnergy))*(dM32*dM32+dm21*dm21+dM32*dm21)
 
    c0 = (-2./27)*A*A*A+(A*A/(18*neutrinoEnergy))*(Ue1*Ue1*(dM32+2*dm21) + Ue2*Ue2*(dM32-dm21) - Ue3*Ue3*(2*dM32+dm21))+(A/(36*neutrinoEnergy*neutrinoEnergy))*(Ue1*Ue1*(2*dM32+dm21)*(dM32-dm21)+Ue2*Ue2*(2*dM32+dm21)*(dM32+2*dm21)-Ue3*Ue3*(dM32+2*dm21)*(dM32-dm21)-(dM32*dM32+dm21*dm21+dM32*dm21))-(1./(216*neutrinoEnergy*neutrinoEnergy*neutrinoEnergy))*(2*dM32+dm21)*(dM32+2*dm21)*(dM32-dm21)
    
    #c0p = -np.linalg.det(T_mass_mat)
    #sg, lnd = np.linalg.slogdet(T_mass_mat)
    
    #Eigenvalues
    #hay un error en las líneas siguientes al utilizar energías del orden de 10 eV, no obstante el cálculo de numpy funciona
    
    s1PlusS2 = 2*np.sqrt((-1./3)*c1)*np.cos((1./3)*np.arctan((1./c0)*np.sqrt(-c0*c0-(4./27)*c1*c1*c1)))
    s1MinusS2 = -2j*np.sqrt((-1./3)*c1)*np.sin((1./3)*np.arctan((1./c0)*np.sqrt(-c0*c0-(4./27)*c1*c1*c1)))
    lam1 = -0.5*s1PlusS2 + np.sqrt(3.)*1j*s1MinusS2/2
    lam2 = -0.5*s1PlusS2 - np.sqrt(3.)*1j*s1MinusS2/2
    lam3 = s1PlusS2
    lam = [lam1, lam2, lam3]
    '''
    lam1 = -np.sqrt((-1./3)*c1)*np.cos((1./3)*np.arctan((1./c0)*np.sqrt(-c0*c0 - (4./27)*c1*c1*c1)))+np.sqrt(-c1)*np.sin(np.arctan((1./c0)*np.sqrt(-c0*c0 - (4./27)*c1*c1*c1)))
    lam2 = -np.sqrt((-1./3)*c1)*np.cos((1./3)*np.arctan((1./c0)*np.sqrt(-c0*c0 - (4./27)*c1*c1*c1)))-np.sqrt(-c1)*np.sin(np.arctan((1./c0)*np.sqrt(-c0*c0 - (4./27)*c1*c1*c1)))
    lam3 = 2*np.sqrt((-1./3)*c1)*np.cos((1./3)*np.arctan((1./c0)*np.sqrt(-c0*c0 - (4./27)*c1*c1*c1)))
    lam = [lam1, lam2, lam3]
    print np.arctan((1./c0)*np.sqrt(-c0*c0 - (4./27)*c1*c1*c1))
'''
    #Calculate operator
    trace_hamiltonian=0.5*E21+E32+3*neutrinoEnergy+A
    phi_phase = np.exp(-1j*L*trace_hamiltonian)
    
    summ=0
    for a in range(3):
        summ+=np.exp(-1j*L*lam[a])*(1./(3*lam[a]*lam[a]+c1))*((lam[a]*lam[a]+c1)*np.identity(3)+lam[a]*T_flav_mat+T_sq_flav_mat)

    summ*=phi_phase
   
    return summ

def calculateProbabilities():
    earth = 1
    sun=0
    if earth:
        coord_init = -6371. #km
        coord_end = 6371. #km
    elif sun:
        coord_init =0. #km
        coord_end = 6.957e5 #km
    
    N =100 #energy steps
    Steps = 100000 #spatial steps
    step_len = np.abs(coord_end-coord_init)/Steps
    
    EnergyLins = np.logspace(1, 13, N)
    
    Probabilities = np.zeros([N,3])
    
    for i in range(N):
        energy = EnergyLins[i]
        coord = coord_init
        operator_product = np.identity(3)
        for k in range(Steps):
            #density = density_to_potential(sun_rho(coord),0)
            density = fig_1_density(coord)
            coord+=step_len
            iter_operator = calculateOperator(energy, density, longitude_units_conversion(step_len))
            operator_product_copy = np.copy(operator_product)
            operator_product = np.matmul(iter_operator, operator_product_copy)
            prob_corr=0
            unit_corr=0
            det_corr=0
            
            if unit_corr:
                unity_check = np.matrix(np.copy(operator_product))
                Id = np.matmul(unity_check, unity_check.H)
                detId = np.linalg.det(Id)
                #print np.abs(detId), '\n\n'
                operator_product/=np.sqrt(detId)
            if prob_corr:
                tot_P=np.abs(operator_product[0,0])**2+np.abs(operator_product[1,0])**2+np.abs(operator_product[2,0])**2
                #print tot_P
                operator_product/=np.sqrt(tot_P)
            if det_corr:
                detU = np.linalg.det(operator_product)
                operator_product/=np.abs(detU)
        for n in range(3):
            Probabilities[i, n] = np.abs(operator_product[n,0])**2
            
    return EnergyLins, Probabilities

#%%
energies, probData = calculateProbabilities()
#%%

#u =calculateOperator(1e4, 1e-13, longitude_units_conversion(100))
#print np.matmul(u, u.H)
#print longitude_units_conversion(100)
#%%
probabilities, ax = plt.subplots(2, 2, figsize=(10, 5))

for i in range(2):
    for k in range(2):
        ax[i,k].set_xscale('log')
        ax[i,k].set_xlabel('$E_{\\nu}$(eV)', fontsize=15)
        ax[i,k].set_xlim(1e3, 1e13)
#ax[1,1].set_ylim(1-0.00001, 1+0.00001)
ax[0,1].set_ylim(0,0.5)
ax[1,0].set_ylim(0,0.5)
ax[0,0].set_ylim(0,1)

ax[0, 0].plot(energies, probData[:,0])
ax[0, 1].plot(energies, probData[:,1])
ax[1, 0].plot(energies, probData[:,2])
ax[1, 1].plot(energies, probData[:,0]+ probData[:,1]+ probData[:,2] )
ax[0, 0].set_ylabel('$P_{e e}$', fontsize=15)
ax[0, 1].set_ylabel('$P_{\mu e}$', fontsize=15)
ax[1, 0].set_ylabel('$P_{\\tau e}$', fontsize=15)
ax[1, 1].set_ylabel('$P_{\\tau e}+P_{\mu e}+ P_{e e}$', fontsize=15)

#ax[1, 1].set_ylim(1-0.01, 1+0.01)

plt.tight_layout()
plt.gcf()
plt.savefig('probPlot.png', dpi=300)
