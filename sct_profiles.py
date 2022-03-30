import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import math

def theta(y):
    if y < 900:
        x = 291
    elif y == 900:
        x = 291
    elif y == 900.01:
        x = 302
    else:
        x = y*(302 -313)/(900 - 2900) + 297

    return x

def qt(y):
    if y < 900:
        x = 10.5
    elif y == 900:
        x = 10.5
    elif y == 900.01:
        x = 4
    else:
        x = y*(3.5 - 4)/(2500 - 900) + 4.5

    return x

def v(y):
    if y < 200:
        x = -5.6
    else:
        x = y*(-1.6 - -5.6)/(2900 - 200) - 5.88

    return x

def pressure_equation(P):
    R = 8.3144598
    g = 9.80665
    M = 0.0289644
    P0 = 101681.1   # Pa
    T0 = 288.15
    L0 = -0.0065
    h0 = 0

   # h = (1/L0) * (T0 * (P/P0)**(-R*L0/g*M) - T0) + h0
    h = (-R*T0)/(M*g)*math.log(P/P0)


    return h

def P_to_h(p_values):
    h_values = np.empty(len(p_values), dtype=float)

    for n,P in enumerate(p_values):
        h_values[n] = pressure_equation(P)

    return h_values

#def sst_increase():

Y = np.asarray([0,450,900,900.01,1000,2000,3000,4000,4250])
Y_v = np.asarray([0,200,200.01,900,1000,2000,3000,4000,4250])
#print(Y)
X_theta = np.empty(len(Y), dtype = float)
X_qt = np.empty(len(Y), dtype = float)
X_v = np.empty(len(Y), dtype = float)

for i in range(len(Y)):
    X_theta[i] = theta(Y[i])
    X_qt[i] = qt(Y[i])
    X_v[i] = v(Y_v[i])
#print(X_qt)

ds = xr.open_dataset('./composite_ref_rrtm_dTs0_dTa0_1xCO2_adjomega.nc')
tak_heights = P_to_h(ds.lev.values)

fig1,ax1 = plt.subplots()
fig2,ax2 = plt.subplots()
fig3,ax3 = plt.subplots()
ax1.plot(X_qt, Y)
ax1.plot(ds.q[0,:,0,0].values*1000, tak_heights)
ax2.plot(X_theta, Y)
ax2.plot(ds.T[0,:,0,0].values, tak_heights)
ax3.plot(X_v, Y_v)
ax1.set_ylim((0,4000))
ax2.set_ylim((0,4000))
ax3.set_xlim((-8.0, 2.0))

ds.Tg.plot()
plt.show()
