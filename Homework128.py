"""HOMEWORK DUE 1/28"""
import numpy as np
import matplotlib.pyplot as plt
mm = 0.001  #mm to m conversion
MPa = 1/(1*(10**6))

'''GIVENS'''
S_ut = 600     #Ultimate tensile strenth MPa
k_a = 0.915
k_b = 0.858
k_c = 0.59
k_d = 1
k_e = 0.897
k_f = 6.76
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
x = (180+20)*mm
M = ((F*(b**2))/(l**3))*(x*(3*a+b)-a*l)

#Stress at point
pi = np.pi
r = D/2
y = r
I = 0.25*pi*r**4
sigma_bend = (M/I)*y

#Total stress
sigma_max = sigma_bend
sigma_min = 0

'''SHEAR'''
#Torque component of shear
J = (pi*D**4)/32
Tau_toruqe_max = (T_max/J)*r
Tau_torque_min = (T_min/J)*r

#Shear at point
R_1 = ((F*b**3)/l**3)*(3*a+b)
V_AB = R_1
Tau_shpt = V_AB

#Total Shear
Tau_max = Tau_toruqe_max+Tau_shpt
Tau_min = Tau_torque_min+Tau_shpt

if Tau_max <= Tau_min:
    print("Max Shear is less than min shear")

"""MEAN AND ALTERNATING STRESSES"""
sigma_mean = (sigma_max+sigma_min)/2
sigma_alt = (sigma_max-sigma_min)/2

tau_mean = (Tau_max+Tau_min)/2
tau_alt = (Tau_max-Tau_min)/2

"""MOHRS CIRCLE"""
#Mean
simga_x = sigma_mean                  #Change this value
simga_y = 0
simga_avg = (simga_x+simga_y)/2
tay_xy = tau_mean                     #Change this value
R = np.sqrt((simga_x-simga_avg)**2+tay_xy**2)

#Mohrs Output
sigma_m1 = simga_avg+R       #Output of the mohrs function
sigma_m2 = simga_avg-R

#Alternating
simga_x = sigma_alt                  #Change this value
simga_y = 0
simga_avg = (simga_x+simga_y)/2
tay_xy = tau_alt                     #Change this value
R = np.sqrt((simga_x-simga_avg)**2+tay_xy**2)

#Mohrs Output
sigma_a1 = simga_avg+R       #Output of the mohrs function
sigma_a2 = simga_avg-R

"""VON MISSES"""
#Mean
Sig1 = sigma_m1      #Change this value
Sig2 = sigma_m2      #Change this value
Sig3 = 0
num = (Sig1-Sig2)**2+(Sig2-Sig3)**2+(Sig3-Sig1)**2
dem = 2
Sig_mean = np.sqrt(num/dem)   #Output

#Alternating
Sig1 = sigma_a1      #Change this value
Sig2 = sigma_a2      #Change this value
Sig3 = 0
num = (Sig1-Sig2)**2+(Sig2-Sig3)**2+(Sig3-Sig1)**2
dem = 2
Sig_alt = np.sqrt(num/dem)   #Output

"""Se"""
S_e_prime = (0.5*S_ut)/MPa
Se = k_a*k_b*k_c*k_d*k_e*k_f*S_e_prime


"""LIFE CYCLES"""

n_inv = (Sig_alt/Se)+(Sig_mean/S_ut)
print(Sig_alt*MPa)
print(Se*MPa)
print(Sig_mean*MPa)
print(Se)
print(1/n_inv)



"""SOLUTION PRINTING"""
print(""*50)
print("PART A SOLUTIONS")
print("Sigma_mean",sigma_mean*MPa,"MPa")
print("Sigma_alternating",sigma_alt*MPa,"MPa")
print("Tau_mean",tau_mean*MPa,"MPa")
print("Tau_alternating",tau_alt*MPa,"MPa")


"""INTERMEDIATE VALUES"""
print(""*50)
print("INTERMEDIATE VALUES")
print("Mohrs output mean (1,2) [Pa]",sigma_m1,sigma_m2)
print("Mohrs output alternating (1,2) [Pa]",sigma_a1,sigma_a2)
print("Principal mean",Sig_mean)
print("Principal alt",Sig_alt)
print("Moment at Point",M,"N-m")
print("a",a)
print("x",x)




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

