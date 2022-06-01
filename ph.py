import numpy
import random
from log import *

print("Choose p, h:")
p, h = map(int, input().split())

GaluaField.set_params(p, h)


a = [0] * GaluaField.p
for i in range(GaluaField.p):
    a[i] = GaluaField.logg(GaluaField.g.to_int(GaluaField.p), Poly([i, 1]).to_int(GaluaField.p))

pi = list(numpy.random.permutation(GaluaField.p))

b = [0] * GaluaField.p
for i in range(GaluaField.p):
    b[i] = a[pi[i]]

d = int(random.random() * (GaluaField.q))

c = [0] * GaluaField.p
for i in range(GaluaField.p):
    c[i] = b[i] + d

print(f"PUBLIC KEY: {c}, {p}, {h}")
print(f"PRIVATE KEY: {GaluaField.f}, {GaluaField.g}, {pi}, {d}")
