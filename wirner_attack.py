import numpy as np


def continued_fraction_algorithm(f_):
    i = 0
    q = rational_to_contfrac(f_)
    while True:
        f = contfrac_to_rational(q[:i + 1])
        if i % 2 == 0:
            tmp = np.zeros(i + 1, dtype=int)
            tmp[i] = 1
            if f == contfrac_to_rational(q[:i + 1] + tmp):
                return f
        else:
            if f == contfrac_to_rational(q[:i + 1]):
                return f
        i += 1



def rational_to_contfrac(f):
    x = f[0]
    y = f[1]
    q = np.array([], dtype=int)
    while y > 0:
        a = x // y
        q = np.append(q, a)
        x, y = y, x - a * y
    return q


def contfrac_to_rational(q):
    m = len(q) - 1
    n0, d0 = 0, 1
    n1, d1 = 1, 0
    for i in range(m + 1):
        n = q[i] * n1 + n0
        d = q[i] * d1 + d0
        n0, d0 = n1, d1
        n1, d1 = n, d
    return [n1, d1]


def wiener_attack(e, N):
    f_ = [e, N]
    f = continued_fraction_algorithm(f_)
    k = f[0]
    d = f[1]
    phi = (e * d - 1) // k
    rem_phi = (e * d - 1) % k
    x = N - phi + 1
    y = (x // 2) ** 2 - N
    if rem_phi == 0 and x % 2 == 0 and is_perfect_square(y):
        return d
    else:
        return None


def is_perfect_square(n):
    n = int(n)
    sq_mod256 = (
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
    if sq_mod256[n & 0xff] == 0:
        return False

    mt = (
        (9, (1, 1, 0, 0, 1, 0, 0, 1, 0)),
        (5, (1, 1, 0, 0, 1)),
        (7, (1, 1, 1, 0, 1, 0, 0)),
        (13, (1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1)),
        (17, (1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1))
    )
    a = n % (9 * 5 * 7 * 13 * 17)
    if any(t[a % m] == 0 for m, t in mt):
        return False

    return isqrt(n) ** 2 == n


def isqrt(n):
    if n == 0:
        return 0

    x = 2 ** ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


if __name__ == '__main__':
    e = 17993
    N = 90581
    d = wiener_attack(e, N)
    print(d)
