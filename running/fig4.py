import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## 稳态对比
data_steady250=pd.read_csv('steady250.csv')
data_250=pd.read_csv('sweep_detuning_250.csv')
fig4=plt.figure()
plt.plot(data_250.delta,data_250.omega,label= 'omega250')
plt.plot(data_250.delta,data_250.n,label= 'n250')
plt.plot(data_250.delta,data_steady250.nbose,label= 'steady_n250')
plt.plot(data_250.delta,data_steady250.pe,label= 'steady_e_femi250')
plt.legend()
plt.show()