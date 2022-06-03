import math
from operator import le
from telnetlib import GA
import numpy
import random
from log import *

print("Choose p, h:")
p, h = map(int, input().split())

GaluaField.set_params(p, h)
#GaluaField.set_params(p, h, Poly([6, 0, 3, 3]), Poly([2, 6, 5, 3, 1]))

a = [0] * GaluaField.p
for i in range(GaluaField.p):
    a[i] = GaluaField.logg(GaluaField.g.to_int(GaluaField.p), Poly([i, 1]).to_int(GaluaField.p))
    print(f'logg ({Poly([i, 1])}) = {a[i]}')
#print(a)

pi = list(numpy.random.permutation(GaluaField.p))

#pi = [6, 4, 0, 2, 1, 5, 3]

b = [0] * GaluaField.p
for i in range(GaluaField.p):
    b[i] = a[pi[i]]

d = int(random.random() * (GaluaField.q))
#d = 1702
c = [0] * GaluaField.p
for i in range(GaluaField.p):
    c[i] = (b[i] + d) % (GaluaField.q - 1)

print(f"PUBLIC KEY: c={c}, p={p}, h={h}")
print(f"PRIVATE KEY: f={GaluaField.f}, g={GaluaField.g}, pi={pi}, d={d}")

pi_ = [0] * GaluaField.p
for i in range(GaluaField.p):
    pi_[pi[i]] = i


print("MESSAGE:")
message = input()
tmp = message

message = ""
for el in tmp:
    binn = bin(ord(el))[2:]
    binn = "0" * (8 - len(binn)) + binn
    message += binn
message_len = len(message)

#print("BINARY VERSION OF MESSAGE IS READY")

C0 = math.factorial(p) // (math.factorial(p - h) * math.factorial(h))
frame_len = int(numpy.log2(C0))


frames_count = message_len // frame_len
if message_len % frame_len != 0:
    frames_count += 1


frames = [0] * frames_count

cur_frames_size = 0
for i in range(0, message_len, frame_len):
    frames[cur_frames_size] = int(message[i : i + frame_len], base=2)
    cur_frames_size += 1

#print("FRAMES ARE READY")

def C(n, k):
    if n < k:
        return 0

    return math.factorial(n) // (math.factorial(n - k) * math.factorial(k))

#print(frames)

m_ready = [None] * frames_count
for k in range(frames_count):
    y = [0] * GaluaField.p
    n = frames[k]
    h_cur = GaluaField.h
    for i in range(1, GaluaField.p + 1):
        C_CUR = C(GaluaField.p - i, h_cur)
        #print(f'C_CUR = {C_CUR}')
        if n >= C_CUR:
            y[i - 1] = 1
            n -= C_CUR
            h_cur -= 1

    m_ready[k] = y

#print("P / H IS READY")
print("MESSAGE ALTERNATIVE FORM:")
print(m_ready)

#c = [12, 11, 6]
encoded = [0] * frames_count
for k in range(frames_count):
    for i in range(GaluaField.p):
        encoded[k] += (m_ready[k][i] * c[i]) % (GaluaField.q - 1)
        encoded[k] %= (GaluaField.q - 1)

print("ENCODED MESSAGE:")
print(encoded)

#decoding

print("DECODING...")

m_dec_ready = [None] * frames_count
for k in range(frames_count):
    print(f'{k}\'th frame...')
    r = (encoded[k] - GaluaField.h * d) % (GaluaField.q - 1)

    print(f'r = ({encoded[k]} - hd) mod {GaluaField.q - 1} = {r}')

    u = GaluaField.pow((GaluaField.g).to_int(GaluaField.p), r)

    print(f'u(x) = g^r = {Poly.to_poly(u, GaluaField.p)}')

    sx = Poly.to_poly(u, GaluaField.p) + GaluaField.f
    sx.normalize_coefficients(GaluaField.p)

    print(f's(x) = u(x) + f = {sx}')
    
    print('s(x) = ', end="")
    y_dec = [0] * GaluaField.p
    for i in range(p):
        value = Poly.to_int(sx, i)
        if value % GaluaField.p == 0:
            root_reverse = (GaluaField.p - i) % GaluaField.p
            print(f'(x + {root_reverse})', end="")
            y_dec[pi_[root_reverse]] = 1
    print()
    m_dec_ready[k] = y_dec

print("DECODED MESSAGE ALTERNATIVE FORM:")
print(m_dec_ready)

decoded = [0] * frames_count
for k in range(frames_count):
    n = 0
    h_cur = GaluaField.h
    for i in range(1, GaluaField.p + 1):
        if m_dec_ready[k][i - 1] == 1:
            n += C(GaluaField.p - i, h_cur)
            h_cur -= 1
    decoded[k] = n

s = ""
cur_len = 0
for i in range(frames_count):
    binn = bin(decoded[i])[2:]
    if i != frames_count - 1 and len(binn) != frame_len:
        binn = "0" * (frame_len - len(binn)) + binn
        cur_len += frame_len
    elif i == frames_count - 1:
        binn = "0" * (message_len - cur_len - len(binn)) + binn
    else:
        cur_len += frame_len
    s += binn
#print(s)

s_s = ((message_len + 7) // 8) * [0]

k = 0
for i in range(0, message_len, 8):
    s_s[k] = chr(int(s[i: i + 8] ,base=2))
    k += 1

print("DECODED MESSAGE:")
print(''.join(s_s)) 