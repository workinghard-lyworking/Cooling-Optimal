import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data=pd.read_csv('eit_sweep_detuning_200.csv')
plt.figure(1)
ax1=plt.subplot(121)
plt.plot(data.delta,data.n,label='n')
plt.legend()
ax2=plt.subplot(122)
plt.plot(data.delta,data.omega_g,label='omega_g')
plt.plot(data.delta,data.omega_r,label='omega_r')
plt.legend()
plt.show()