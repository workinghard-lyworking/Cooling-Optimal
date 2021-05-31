import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 行波边带冷却
## 边带与最优失谐对比

evo1=pd.read_csv('250_-1.CSV')
evo2=pd.read_csv('250_-0.845.CSV')
fig1=plt.figure()
plt.plot(evo1.t,evo1.evo)
plt.plot(evo2.t,evo2.evo)
plt.show()
