import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## 稳态对比
data_steady240=pd.read_csv('steady240.csv')
data_240=pd.read_csv('sweep_detuning_240.csv')
fig4=plt.figure()
plt.plot(data_240.delta,data_240.omega,label= 'omega240')
plt.plot(data_240.delta,data_240.n,label= 'n240')
plt.plot(data_240.delta,data_steady240.nbose,label= 'steady_n240')
plt.plot(data_240.delta,data_steady240.pe,label= 'steady_e_femi240')
plt.legend()
plt.show()