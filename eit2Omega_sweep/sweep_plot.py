import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data3d=pd.read_csv('omegag_omegar_sweep200_65.csv')
nsurf=np.array(data3d.iloc[:,2:])
fig3=plt.figure()
#ax = fig3.gca(projection='3d')
X, Y = np.meshgrid(data3d.omega_g, data3d.omega_r)
#ax.plot_surface(X,Y,nsurf)
plt.contour(X,Y,nsurf,20)
plt.show()
fig3.savefig('test.eps',dpi=600,format='eps')
