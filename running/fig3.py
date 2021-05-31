import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# # 3d+最优曲线
data3d=pd.read_csv('delta_omega_sweep.csv')
nesurf=np.array(data3d.iloc[:,2:])
row=int(nesurf.shape[1]/2)
nsurf=nesurf[:,0:row]
esurf=nesurf[:,row:]
fig3=plt.figure()
ax = fig3.gca(projection='3d')
X, Y = np.meshgrid(data3d.omega, data3d.delta)
ax.plot_surface(X,Y,nsurf)
plt.show()