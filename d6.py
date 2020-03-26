import numpy as np
c = 299792.458  # Ljoshradi km/s
r0 = [0, 0, 6370, 0.0001]  # hnit mottakara
rho = 26570  # km

# Aðferð til að fá hnit dreifð jafnt yfir jörðina
def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec.T

R = sample_spherical(8)
