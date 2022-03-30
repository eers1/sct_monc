import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts

def find_volume(N_mode, r, sigma):
    V = (4*np.pi/3)*N_mode*r**3*np.exp(9*(np.log(sigma)**2/2))
    return V

def find_mass(N_mode, r, sigma, rho):
    V = find_volume(N_mode, r, sigma)
    m = V*rho
    return m

def dist(N0, r, sigma, mu):
    n_l = np.empty(len(r))
    for i,er in enumerate(r):
        n_l[i] = (N0/((2*np.pi)**0.5)*sigma)*np.exp(-(np.log(er)-mu)**2/2*sigma**2)
    return n_l

def create_list(len_bl, len_ft, bl, ft):
    bl=[bl]*len_bl
    ft=[ft]*len_ft
    full=bl+ft
    return full

#m_low = find_mass(100e6, 0.1e-6,1.5, 1800)
#m_high = find_mass(200e6, 0.2e-6,1.5,2600)
#m_mid = find_mass(150e6, 0.15e-6,1.5, 2000)

#m_aitken = find_mass(300e6, 0.1e-7, 2000)

#  schulze et al gives a table of values used as model constraints D
ait_N_bl=150e6
ait_N_ft=200e6
ait_r=0.5*0.05e-6
ait_s=1.25
#ait_rho=1500

acc_N_bl=500e6 #150e6
acc_N_ft=500e6 #100e6
acc_r=0.5*0.2e-6
acc_s=1.5
rho=1500  #kg m^-3

aitken_bl = find_mass(ait_N_bl, ait_r, ait_s, rho)
aitken_ft = find_mass(ait_N_ft, ait_r, ait_s, rho)

accum_bl = find_mass(acc_N_bl, acc_r, acc_s, rho)
accum_ft = find_mass(acc_N_ft, acc_r, acc_s, rho)

accum_source = find_mass(70e6, 0.2e-6, 1.5, 1500)

print("{:.2e}".format(accum_bl), accum_ft)
print(aitken_bl, aitken_ft)
print(accum_source)

len_bl=2 #19
len_ft=42
ait_N = create_list(len_bl,len_ft,ait_N_bl,ait_N_ft)
acc_N = create_list(len_bl,len_ft,acc_N_bl,acc_N_ft)
ait_m = create_list(len_bl,len_ft,aitken_bl,aitken_ft)
acc_m = create_list(len_bl,len_ft,accum_bl,accum_ft)
# ait_N = create_list(len_bl,len_ft,"{:.2e}".format(ait_N_bl),"{:.2e}".format(ait_N_ft))
# acc_N = create_list(len_bl,len_ft,"{:.2e}".format(acc_N_bl),"{:.2e}".format(acc_N_ft))
# ait_m = create_list(len_bl,len_ft,"{:.2e}".format(aitken_bl),"{:.2e}".format(aitken_ft))
# acc_m = create_list(len_bl,len_ft,"{:.2e}".format(accum_bl),"{:.2e}".format(accum_ft))

print(ait_N)
print(len(acc_N), acc_N)
print(ait_m)
print(len(acc_m), acc_m)
#ri=np.linspace(2e-9,10e-6,100000)
ri=np.linspace(1e-9,1000e-9,100000)

fig, ax = plt.subplots()
ait_bl=dist(ait_N_bl, ri, ait_s, np.log(ait_r))
acc_bl=dist(acc_N_bl, ri, acc_s, np.log(acc_r))
bl=dist(acc_N_bl+ait_N_bl, ri, 1.35, np.log(ait_r))
bl = ait_bl+acc_bl
ait_ft=dist(ait_N_ft, ri, ait_s, np.log(ait_r))
acc_ft=dist(acc_N_ft, ri, acc_s, np.log(acc_r))
ft = ait_ft+acc_ft
ait_t=dist(ait_N_bl+ait_N_ft, ri, ait_s, np.log(ait_r))
acc_t=dist(acc_N_bl+acc_N_ft, ri, acc_s, np.log(acc_r))
t = ait_t+acc_t
ax.plot(ri*1e9, ait_bl,c='C0') # r in nm, dN in cm^-3
ax.plot(ri*1e9, acc_bl,c='C0', label='Boundary layer')
ax.plot(ri*1e9,bl,'--',c='C0',alpha=0.5)#, label='Boundary layer')
ax.plot(ri*1e9, ait_ft,c='C1') # r in nm, dN in cm^-3
ax.plot(ri*1e9, acc_ft,c='C1', label='Free troposphere')
ax.plot(ri*1e9,ft,'--',c='C1',alpha=0.5)#, label='Free troposphere')
ax.plot(ri*1e9, ait_t,c='C2') # r in nm, dN in cm^-3
ax.plot(ri*1e9, acc_t,c='C2', label='Total')
ax.plot(ri*1e9,t,'--',c='C2',alpha=0.5)#, label='Total')
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('r (nm)')
ax.set_ylabel('dN/dlog(r)')
#ax.set_ylim(1e-1,300)
plt.legend()
plt.show()

