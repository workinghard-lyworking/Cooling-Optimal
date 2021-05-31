import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data=pd.read_csv('eit_four_opt.csv')
plt.figure()
plt.plot(data.list_T,data.list_ZERO,label='initial')
plt.plot(data.list_T,data.list_40,label='t=40')
plt.plot(data.list_T,data.list_90,label='t=90')
plt.plot(data.list_T,data.list_150,label='t=150')
plt.plot(data.list_T,data.list_250,label='t=250')
plt.legend()
plt.show()