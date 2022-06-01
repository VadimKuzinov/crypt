from log import *

print("Choose p, h:")
p, h = map(int, input().split())

GaluaField.set_params(p, h)
print(GaluaField.f)
print(GaluaField.g)

print('logg = ', GaluaField.logg(GaluaField.g.to_int(GaluaField.p), 7))