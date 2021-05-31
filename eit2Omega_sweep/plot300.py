import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data=pd.read_csv('eit_sweep_detuning_300.csv')
datare=pd.read_csv('eit_sweep_detuning_300re.csv')
plt.figure(1)
ax1=plt.subplot(121)
plt.plot(data.delta,data.n,label='n')
plt.plot(datare.delta,datare.n,label='n_re')
plt.legend()
ax2=plt.subplot(122)
plt.plot(data.delta,data.omega_g,label='omega_g')
plt.plot(data.delta,data.omega_r,label='omega_r')
plt.plot(datare.delta,datare.omega_g,label='omega_g_re')
plt.plot(datare.delta,datare.omega_r,label='omega_r_re')
plt.legend()
plt.show()
