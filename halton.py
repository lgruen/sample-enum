#!/usr/bin/env python

primes = [] # prime numbers
faure, faure_inverse = [], [] # Faure permutations

def is_prime(n):
    """Very simple deterministic primality check."""

    for i in range(2, n):
        if not n % i:
            return False

    return True

def init_primes(num_primes):
    """Generate num_primes prime numbers."""

    global primes
    primes = []
    candidate = 2
    while num_primes:
        if is_prime(candidate):
            primes.append(candidate)
            num_primes -= 1
        candidate += 1

def get_faure_permutation(b):
    """Compute the Faure digit permutation for 0, ..., b - 1."""

    if b < 2:
        return (0,)

    elif b == 2:
        return (0, 1)

    elif b & 1: # odd
        c = (b - 1) / 2

        def faure_odd(i):
            if i == c:
                return c

            f = faure[b - 1][i - int(i > c)]
            return f + int(f >= c)

        return tuple((faure_odd(i) for i in range(b)))

    else: # even
        c = b / 2

        def faure_even(i):
            if i < c:
                return 2 * faure[c][i]
            else:
                return 2 * faure[c][i - c] + 1

        return tuple((faure_even(i) for i in range(b)))

def invert_permutation(perm):
    """Returns the inverted permutation of 0, ..., len(perm) - 1."""

    inverse = [0] * len(perm)
    for (i, j) in enumerate(perm):
        inverse[j] = i

    return tuple(inverse)

def init_faure_permutations():
    """Initialize the list of Faure permutations up to the maximum prime number."""

    global faure, faure_inverse
    faure, faure_inverse = [], []
    for b in range(primes[-1] + 1):
        faure.append(get_faure_permutation(b))
        faure_inverse.append(invert_permutation(faure[b]))

def get_multiplicative_inverses(n):
    """Initialize the list of multiplicative inverses, by
    computing ((N / n[i])^(-1) mod n[i]) for
    N = n[0] * n[1] * ... * n[len(n) - 1]."""
    
    n_product = product(n)

    return tuple((multiplicative_inverse(n_product // n_i, n_i) for n_i in n))

def init(num_components):
    """Initializes the precomputed lists to allow for up to
    num_components-dimensional samples."""

    init_primes(num_components)
    init_faure_permutations()

def halton(index, base):
    """Returns the Halton radical inverse of the integer index
    in the given base. The result is a tuple (inverse, num_digits),
    where inverse is an integer consisting of num_digits flipped
    digits of index. The floating point representation of the
    radical inverse is given by (inverse / base ** num_digits)."""

    assert isinstance(index, (int, long)) and index >= 0
    assert isinstance(base, (int, long)) and base >= 2

    (index, digit) = divmod(index, base)
    inverse = faure[base][digit]
    num_digits = 1
    while index:
        (index, digit) = divmod(index, base)
        inverse = inverse * base + faure[base][digit]
        num_digits += 1

    return (inverse, num_digits)

def halton_inverse((inverse, num_digits), base):
    """Given a result (inverse, num_digits) from the function
    halton, this function reconstructs the index. In other
    words, halton_inverse(halton(index, base), base) == index."""

    assert isinstance(inverse, (int, long)) and inverse >= 0
    assert isinstance(num_digits, (int, long)) and num_digits >= 1
    assert isinstance(base, (int, long)) and base >= 2

    index = 0
    while num_digits:
        (inverse, digit) = divmod(inverse, base)
        index = index * base + faure_inverse[base][digit]
        num_digits -= 1

    return index

def extended_euclid(a, b):
    """Returns a triple of the form (d, x, y) that satisfies
    d = gcd(a, b) = ax + by, where gcd(a, b) is the greatest
    common divisor of a and b. Note that x and y may be negative."""

    assert isinstance(a, (int, long)) and a >= 0
    assert isinstance(b, (int, long)) and b >= 0

    if b == 0:
        return (a, 1, 0)

    (div, mod) = divmod(a, b)
    (d, x, y) = extended_euclid(b, mod)

    return (d, y, x - div * y)

def multiplicative_inverse(a, n):
    """Compute the multiplicative inverse of a, modulo n,
    when a and n are relatively prime. The result is in the
    range 1, ..., n - 1."""

    assert isinstance(a, (int, long)) and a > 0
    assert isinstance(n, (int, long)) and n > 1

    (d, x, y) = extended_euclid(a, n)
    return x % n

def product(n):
    """Return n[0] * n[1] * ... * n[len(n) - 1]."""
    
    product = 1
    for p in n:
        product *= p

    return product

def chinese_remainder(a, n, multiplicative_inverses):
    """a and n are lists of integers of the same length k.
    If n[0], n[1], ..., n[k - 1] are pairwise relatively prime,
    then for any integers a[0], a[1], ..., a[k - 1] the set
    of simultaneous equations x = a[i] % n[i], for i = 0, ..., k - 1,
    has a unique solution modulo N = n[0] * n[1] * ... * n[k - 1].
    multiplicative_inverses must be a list of (m_i^(-1) mod n[i]),
    where m_i = N / n[i].
    This function returns the unknown x, based on the Chinese
    remainder theorem."""

    assert len(a) == len(n) and len(n) > 0

    n_product = product(n)

    x = 0
    for (i, a_i) in enumerate(a):
        m_i = n_product // n[i]
        x += a_i * m_i * multiplicative_inverses[i]

    return x % n_product

def get_prime_powers(s):
    """Given a list of integers s[0], s[1], ..., s[k],
    return a tuple of prime powers (p_0 ** m_0, p_1 ** m_1, ..., p_k ** m_k)
    such that p_0 ** m_0 >= s[0], p_1 ** m_1 >= s[1], ..., p_k ** m_k >= s[k],
    together with the exponents (m_0, m_1, ..., m_k)."""
    
    assert len(s) <= len(primes)

    prime_powers, exponents = [], []
    for (i, n) in enumerate(s):
        value, exponent = primes[i], 1
        while value < n:
            value *= primes[i]
            exponent += 1

        prime_powers.append(value)
        exponents.append(exponent)

    return (tuple(prime_powers), tuple(exponents))

def get_component(index, k):
    """Return the k-th component of the Halton-sequence sample for the given
    index. The result is a floating-point number in [0, 1)."""

    s = halton(index, primes[k])
    return float(s[0]) / primes[k] ** s[1]

def get_sample(index, num_components):
    """Return the num_components-dimensional Halton sample for the given index.
    The result is in the num_components-dimensional unit-cube."""

    return tuple((get_component(index, i) for i in range(num_components)))

def get_scaled_component(index, k, prime_exponent):
    """Return the k-th component of the Halton-sequence sample for the given
    index. The result is scaled to the voxel grid.."""

    s = halton(index, primes[k])
    return float(s[0]) / primes[k] ** (s[1] - prime_exponent)

def get_scaled_sample(index, num_components, prime_exponents):
    """Return the num_components-dimensional Halton sample for the given index,
    scaled by the prime powers. That means the result lies inside the voxel grid."""
    
    return tuple((get_scaled_component(index, i, prime_exponents[i])
                 for i in range(num_components)))

def get_offset(coords, (prime_powers, prime_exponents), multiplicative_inverses):
    """Return the fixed offset of the index for the Halton sequence for a
    voxel stratum identified by the given coordinates."""

    assert len(coords) == len(prime_exponents)

    inverted = tuple((halton_inverse((c, prime_exponents[i]), primes[i])
               for (i, c) in enumerate(coords)))

    return chinese_remainder(inverted, prime_powers, multiplicative_inverses)

def test_sampling():
    """Test the functionality by generating samples in a 3D-voxel grid."""

    res = (8, 10, 9) # voxel resolution
    (prime_powers, prime_exponents) = get_prime_powers(res)
    multiplicative_inverses = get_multiplicative_inverses(prime_powers)

    increment = product(prime_powers)

    num_samples_per_voxel = 4
    num_components = 3

    # iterate over voxels
    for x in range(res[0]):
        for y in range(res[1]):
            for z in range(res[2]):
                offset = get_offset((x, y, z),
                                    (prime_powers, prime_exponents),
                                    multiplicative_inverses)

                for i in range(num_samples_per_voxel):
                    index = offset + i * increment
                    s = get_scaled_sample(index, num_components, prime_exponents)

                    assert x <= s[0] < x + 1
                    assert y <= s[1] < y + 1
                    assert z <= s[2] < z + 1
        
# Run some tests if executed as stand-alone program.
if __name__ == "__main__":
    max_num_components = 10
    init(max_num_components)

    assert faure[7] == (0, 2, 5, 3, 1, 4, 6)
    assert faure[8] == (0, 4, 2, 6, 1, 5, 3, 7)
    assert faure_inverse[7] == (0, 4, 1, 3, 5, 2, 6)
    assert faure_inverse[8] == (0, 4, 2, 6, 1, 5, 3, 7)

    for i in range(2 ** 4):
        assert halton_inverse(halton(i, 2), 2) == i

    for i in range(3 ** 3):
        assert halton_inverse(halton(i, 3), 3) == i

    for i in range(1, 31):
        j = multiplicative_inverse(i, 31)
        assert i * j % 31 == 1

    a = (123, 456, 7891)
    n = (2 ** 10, 3 ** 7, 5 ** 8)
    multiplicative_inverses = get_multiplicative_inverses(n)
    x = chinese_remainder(a, n, multiplicative_inverses)
    for (i, a_i) in enumerate(a):
        assert x % n[i] == a[i]

    res = (640, 480, 512)
    (prime_powers, prime_exponents) = get_prime_powers(res)
    for (i, k) in enumerate(res):
        assert prime_powers[i] >= k
        assert primes[i] ** prime_exponents[i] == prime_powers[i]

    test_sampling()

