import numpy as np
gamma=20
nu=1.0
def nss(Delta,Omegar):
    Ap=1/(gamma**2*nu**2+4*(Omegar**2/4-nu*(nu-Delta))**2)
    Am=1/(gamma**2*nu**2+4*(Omegar**2/4-nu*(nu+Delta))**2)
    return Ap/(Am-Ap)

print(nss(74.6562,np.sqrt(2)*np.sqrt(2+np.sqrt(gamma**2+4*74.6562**2))))