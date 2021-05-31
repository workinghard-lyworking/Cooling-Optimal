import matplotlib.pyplot as plt
import numpy as np

gamma=0.1
alpha=0.4
nu=1.0
def nss(Delta):
    Ap=gamma/((Delta-nu)**2+(gamma**2)/4)+alpha*gamma/(Delta**2+(gamma**2)/4)
    Am=gamma/((Delta+nu)**2+(gamma**2)/4)+alpha*gamma/(Delta**2+(gamma**2)/4)
    return Ap/(Am-Ap)

x=np.arange(-1.1,-0.6,0.01)
fig=plt.figure()
plt.plot(x,nss(x))
print(nss(-0.8905))
plt.show()
