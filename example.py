import cplxfnc

s=1+1j
a=3-5j
zeta = cplxfnc.zeta(s, a)
print("zeta({}, {}) = {}".format(s, a, zeta))

s=-0.1
z=-3.4
gamma_inc = cplxfnc.gamma_inc(s, z)
print("gamma_inc({}, {}) = {}".format(s, z, gamma_inc))

gamma_inc = cplxfnc.gamma_inc(s, 0)
