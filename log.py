import numpy

class Poly:
    x: list = []
    deg: int = 0

    def __init__(self, x: list, deg: int = None):
        self.x = x.copy()
        self.deg = deg or len(x)
    
    def __str__(self):
        return ' + '.join(['{0}x^{1}'.format(self.x[i], i) for i in range(self.deg) if self.x[i] != 0])
    
    def __mul__(self, other):
        result = Poly([0] * (self.deg + other.deg - 1), self.deg + other.deg - 1)
        for i in range(self.deg):
            for j in range(other.deg):
                result.x[i + j] += self.x[i] * other.x[j]
        return result

    def __add__(self, other):
        if self.deg < other.deg:
            self, other = other, self

        result = Poly(self.x, self.deg)
        for i in range(other.deg):
            result.x[i] += other.x[i]

        for i in range(self.deg - 1, -1, -1):
            if result.x[i] == 0:
                del result.x[i]
                result.deg -= 1
            else:
                break

        return result

    def __sub__(self, other):
        is_neg = False
        if self.deg < other.deg:
            self, other = other, self
            is_neg = True

        result = Poly(self.x, self.deg)
        for i in range(other.deg):
            result.x[i] -= other.x[i]
            
        if is_neg:
            for i in range(result.deg):
                result.x[i] *= (-1)

        for i in range(self.deg - 1, -1, -1):
            if result.x[i] == 0:
                del result.x[i]
                result.deg -= 1
            else:
                break

        return result

    def __mod__(self, other):
        current_poly = Poly(self.x, self.deg)
        if other is None:
            return current_poly

        #current_poly = Poly(self.x, self.deg)
        while current_poly.deg >= other.deg:

            mul_deg = current_poly.deg - other.deg            
            mul_kf = current_poly.x[-1] / other.x[-1]
            mul_poly = Poly([0] * (mul_deg + 1), mul_deg + 1)
            mul_poly.x[-1] = mul_kf
            dif_poly = other * mul_poly
            current_poly = current_poly - dif_poly
        return current_poly

    def __floordiv__(self, other):
        if other is None:
            return self

        result = Poly([])
        current_poly = Poly(self.x, self.deg)
        while current_poly.deg >= other.deg:
            mul_deg = current_poly.deg - other.deg
            mul_kf = current_poly.x[-1] / other.x[-1]
            mul_poly = Poly([0] * (mul_deg + 1), mul_deg + 1)
            mul_poly.x[-1] = mul_kf
            result = result + mul_poly
            dif_poly = other * mul_poly
            current_poly = current_poly - dif_poly

        return result

    def __pow__(self, exp: int):
        result = Poly([1])
        if exp == 0:
            return result

        for i in range(exp):
            result = result * self
        return result

    def __eq__(self, other):
        flag = True
        if self.deg != other.deg:
            flag = False

        if not flag:
            return False

        for i in range(self.deg):
            if self.x[i] != other.x[i]:
                flag = False
                break
        return flag

    def normalize_coefficients(self, mod: int):
        
        for i in range(self.deg):
            self.x[i] %= mod

        old_deg = self.deg
        for i in range(old_deg - 1, -1 ,-1):
            if self.x[i] == 0:
                self.deg -= 1
                del self.x[i]
            else:
                break

    def to_int(self, base: int):
        res = 0
        exp = 1
        for i in range(self.deg):
            res += exp * self.x[i]
            #print(res)
            exp *= base
        return res

    @staticmethod
    def to_poly(num: int, base: int):
        result = Poly([], 0)
        while num > 0:
            digit = num % base
            num //= base
            result.x.append(digit)
            result.deg += 1

        return result


class GaluaField:
    p: int = 0
    h: int = 0
    q: int = 0
    f: Poly = None
    g: Poly = None
    factorization: list = None

    @classmethod
    def set_params(cls, p, h, g=None, f=None):
        cls.p = p
        cls.h = h
        cls.g = g
        cls.f = f
        cls.q = p**h
        cls.factorization = get_factorization(cls.q - 1)
        if f is None:
            cls.f = next(cls.get_irreducible_poly())
        if g is None:
            cls.g = next(cls.get_primitive_poly())

    @classmethod
    def set_irreducible_poly(cls, poly: Poly):
        cls.f = poly

    @classmethod
    def set_primitive_poly(cls, poly: Poly):
        cls.g = poly

    @classmethod
    def multiply(cls, a: int, b: int):
        result = Poly.to_poly(a, cls.p) * Poly.to_poly(b, cls.p)
        result.normalize_coefficients(GaluaField.p)
        return cls.mod(result.to_int(cls.p), cls.f.to_int(cls.p))
        result = result % cls.f
        #result = cls.mod(result.to_int(cls.p), cls.f.to_int(cls.p))
        result.normalize_coefficients(cls.p)

        #return result
        return result.to_int(cls.p)

    @classmethod
    def divmod(cls, a: int, b: int):
        result = Poly([])
        current_poly = Poly.to_poly(a, cls.p)
        b = Poly.to_poly(b, cls.p)
        b.normalize_coefficients(cls.p)
        while current_poly.deg >= b.deg:
            mul_deg = current_poly.deg - b.deg
            while (current_poly.x[-1] % b.x[-1]) != 0:
                current_poly.x[-1] += cls.p

            mul_kf = current_poly.x[-1] // b.x[-1]
            mul_poly = Poly([0] * (mul_deg + 1), mul_deg + 1)
            mul_poly.x[-1] = mul_kf
            result = result + mul_poly
            dif_poly = b * mul_poly
            current_poly = current_poly - dif_poly

        #result = result % cls.f
        result.normalize_coefficients(cls.p)

        current_poly.normalize_coefficients(cls.p)

        #print(f'{Poly.to_poly(a, cls.p)} = ({result}) * ({b}) + {current_poly}')
        return (result.to_int(cls.p), current_poly.to_int(cls.p))

    @classmethod
    def div(cls, a: int, b: int):
        return cls.divmod(a, b)[0]

    @classmethod
    def mod(cls, a: int, b: int):
        return cls.divmod(a, b)[1]

    @classmethod
    def add(cls, a: int, b: int):
        result = Poly.to_poly(a, cls.p) + Poly.to_poly(b, cls.p)
        result = result % cls.f
        result.normalize_coefficients(cls.p)
        return result.to_int(cls.p)

    @classmethod
    def sub(cls, a: int, b: int):
        result = Poly.to_poly(a, cls.p) - Poly.to_poly(b, cls.p)
        result = result % cls.f
        result.normalize_coefficients(cls.p)
        return result.to_int(cls.p)

    @classmethod
    def get_inverse(cls, a: int):
        t = 0
        newt = 1
        r = cls.f.to_int(cls.p) 
        newr = a

        while newr != 0:
            quotient = cls.div(r, newr)     
            r, newr = newr, cls.sub(r, cls.multiply(quotient, newr))
            t, newt = newt, cls.sub(t, cls.multiply(quotient, newt))

        if r >= cls.p:
            return -1

        k = 0
        while (1 + cls.p * k) % r != 0:
            k += 1
        
        newk = 1 + cls.p * k
        newk = newk // r
        result = cls.multiply(newk, t)
        return result

    @classmethod
    def pow(cls, x: int, exp: int):
        #print(f'{Poly.to_poly(x, cls.p)} ^ {exp} by mod {cls.f}')
        if exp < 0:
            x = cls.get_inverse(x)
            exp = -exp
        t = 1
        for i in range(exp):
            t = cls.multiply(t, x)
            #print(Poly.to_poly(t, cls.p))
        return t

    @classmethod
    def CTOR(cls, xi: list, m: list):
        M = 1
        for el in m:
            M *= el

        x0 = 0
        for i in range(len(xi)):
            MI = M // m[i]
            MI_ = get_inverse_by_mod(MI, m[i])
            x0 = (x0 + xi[i] * MI * MI_) % M
        return x0

    @classmethod
    def logg(cls, base: int, value: int):
        a = base
        b = value
        xi = []
        for pair in cls.factorization:
            y = []

            v = cls.pow(a, (cls.q - 1) // pair[0])
            exp_v = []
            exp = 1
            for i in range(pair[0]):
                exp_v.append(exp)
                exp = cls.multiply(exp, v)
            #print(exp_v)

            power = (cls.q - 1) // pair[0]
            b_exp = cls.pow(b, power)
            #print(b_exp)
            flag = False
            for i in range(pair[0]):
                if exp_v[i] == b_exp:
                    y.append(i)
                    flag = True
                    break
            if not flag:
                return -1

            power_a = 0
            p_exp = 1
            for i in range(1, pair[1]):
                power //= pair[0]
                power_a -= y[i - 1] * p_exp
                p_exp *= pair[0]
                
                a_exp = cls.pow(a, power_a)
                ba = cls.multiply(b, a_exp)
                ba_exp = cls.pow(ba, power)
                #print(ba_exp)
                flag = False
                for k in range(pair[0]):
                    if exp_v[k] == ba_exp:
                        y.append(k)
                        flag = True
                        break
                if not flag:
                    return -1
            
            xi_cur = 0
            p_exp = 1
            for i in range(pair[1]):
                xi_cur += p_exp * y[i]
                p_exp *= pair[0]
            xi.append(xi_cur)
            #print(xi)

        return cls.CTOR(xi, [pair[0] ** pair[1] for pair in cls.factorization])

    @classmethod
    def gcd(cls, x: int, y: int):
        #print(f'gcd({str(Poly.to_poly(x, cls.p))}, {str(Poly.to_poly(y, cls.p))})')
        if x < y:
            x, y = y, x

        #print(f'x={x} y={y}')
        if y == 0:
            return x

        return cls.gcd(y, cls.mod(x, y))

    @classmethod
    def get_irreducible_poly(cls, deg: int = None):
        tmp = cls.f

        p = cls.p
        m = deg or cls.h

        exp = p ** m

        for cur in numpy.random.permutation(exp):
            if cur % cls.p == 0:
                continue

            flag = False 
            cur_poly = exp + cur
            #print("CURRENT_POLY: ", Poly.to_poly(cur_poly, cls.p))
            u = p
            for i in range(m // 2):
                cls.f = Poly.to_poly(cur_poly, cls.p)
                u = cls.pow(u, p)
                cls.f = None
                d = cls.gcd(cur_poly, cls.sub(u, p))
                #print("GCD=", d)
                if d >= cls.p:
                    flag = True
                    break
            if not flag:
                cls.f = tmp
                #print("SUCCESS")
                yield Poly.to_poly(cur_poly, cls.p)

    @classmethod
    def get_irreducible_poly_list(cls, deg: int = None):
        res = []
        tmp = cls.f

        p = cls.p
        m = deg or cls.h

        exp = p ** m
        cur = 0
        flag = True
        while cur < exp - 1:
            flag = False 
            cur += 1
            cur_poly = exp + cur
            u = p
            #print(Poly.to_poly(cur_poly, cls.p))
            for i in range(m // 2):
                cls.f = Poly.to_poly(cur_poly, cls.p)
                u = cls.pow(u, p)
                cls.f = None
                d = cls.gcd(cur_poly, cls.sub(u, p))
                if d >= cls.p:
                    flag = True
                    break
            cls.f = tmp
            if not flag:
                #cls.f = tmp
                res.append(Poly.to_poly(cur_poly, cls.p))
        return res

    @classmethod
    def get_primitive_poly(cls, deg: int = None):
        for degree in numpy.random.permutation(cls.h - 1):
            m  = degree + 1
            flag = True
            for newf in cls.get_irreducible_poly(m):
                #print("NEWF = ", newf)
                flag = False
                for pair in cls.factorization:
                    exp = (cls.q - 1) // pair[0]
                    l = cls.pow(newf.to_int(cls.p), exp)
                    #print(f'{Poly.to_poly(cls.p, cls.p)} ^ {exp} = {Poly.to_poly(l, cls.p)} MOD {newf}')
                    if l == 1:
                        flag = True
                        break
                if not flag:
                    yield newf

def get_factorization(num: int):
    result = list()
    p = num
    div = 2
    cur_exp = 0
    while p != 1:
        if p % div == 0:
            p //= div
            cur_exp += 1
        elif p == 1 or cur_exp != 0:
            result.append((div, cur_exp))
            div += 1
            cur_exp = 0
        else:
            div += 1
    result.append((div, cur_exp))
    return result

def get_inverse_by_mod(x: int, mod: int):
    r = mod
    newr = x
    t = 0
    newt = 1
    
    while newr != 0:
        quotient = r // newr
        t, newt = newt, t - quotient * newt
        r, newr = newr, r - quotient * newr
    
    if r > 1:
        return -1

    if t < 0:
        t = t + mod

    return t

