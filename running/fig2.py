import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## 多个终点时间演化

evo100=pd.read_csv('sweep_detuning_100.csv')
evo250=pd.read_csv('sweep_detuning_250.csv')
evo400=pd.read_csv('sweep_detuning_400.csv')
fig2=plt.figure()
plt.plot(evo100.delta,evo100.n);plt.plot(evo100.delta,evo100.omega)
plt.plot(evo250.delta,evo250.n);plt.plot(evo250.delta,evo250.omega)
plt.plot(evo400.delta,evo400.n);plt.plot(evo400.delta,evo400.omega)
plt.show()