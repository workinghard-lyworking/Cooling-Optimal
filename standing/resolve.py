
gamma=0.1
nu=1.0
def nss(Delta):
    Ap=gamma/((Delta-nu)**2+(gamma**2)/4)
    Am=gamma/((Delta+nu)**2+(gamma**2)/4)
    return Ap/(Am-Ap)

print(nss(-1))