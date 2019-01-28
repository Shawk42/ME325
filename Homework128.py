"""HOMEWORK DUE 1/28"""
import numpy as np
import matplotlib.pyplot as plt
mm = 0.001  #mm to m conversion


'''GIVENS'''
S_ut = 600     #Ultimate tensile strenth MPa
k_a = 0.915
k_b = 0.858
k_c = 0.59
k_d = 1
k_e = 0.897
k_f = 676
T_max = 150 #N-m Max torque exterted on shaft
T_min = 50  #N-m Min torque exterted on shaft
F = 6000    #N
D = 35*mm   #Diamter of shaft at point of interest [m]

'''BENDING STRESS'''
#To calculate the bending stress a shear and moment diagram is needed to find the moment
l = 500*mm   #Length of shaft in meters
b = 175*mm   #Distance from right end of shaft to force [m]
a = l-b      #Distance from left end of shaft to force [m]

#Moment at Point
x = 180*mm
M = ((F*b**2)/l**3)*(x*(3*a+b)-a*l)

#Moment Diagram Plotting
x_ab = np.linspace(0,a,10000)
x_bc = np.linspace(a,l,10000)
M_AB = ((F*b**2)/l**3)*(x_ab*(3*a+b)-a*l)
M_BC = M_AB-F*(x_bc-a)
M_all = np.append(M_AB,M_BC)
x_beam = np.append(x_ab,x_bc)
plt.plot(x_beam,M_all)
plt.plot(x,M,"*")
plt.xlabel("Meters")
plt.ylabel("Moment (N-m)")
plt.grid()
plt.show()

#Stress at point
pi = np.pi
r = D/2
y = r
I = 0.25*pi*r**4

sigma_bend = (M/I)*y




