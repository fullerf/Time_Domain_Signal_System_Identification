from __future__ import division, print_function
from numpy import (zeros, zeros_like, ones_like, 
    complex128, float64, dot, argwhere, where, sqrt, log, pi, abs,
    arange, exp, linspace)
from scipy.linalg.lapack import zggev
import itertools as it

__all__ = ['generate_u', 'gen_amplitudes', 'solve_fdm', 'make_lorenzian', 
        'make_basis'] 

try:
    from _fdm import accurate_pow
except ImportError:
    def accurate_pow(z, n):
        '''perform fewer multiplications to calculate an integer power
        '''
        assert int(n) == n, 'only integer powers supported'
        if n < 0:
            return 1.0/accurate_pow(z, -n)
        else:
            res = 1;
            while n > 1:
                if n % 2 == 1:
                    res *= z
                z *= z
                n = n >> 1 # bitwise shift
            if n > 0:
                res *= z
            return res

def generate_u(ul, signal, gen_u2=False): # generate G0, G1
    '''generate $U^0$ and $U^1$ matricies required for the method
    '''
    NPOW = 3
    assert (len(ul.shape) == 1 and len(signal.shape) == 1), "expected 1d vectors"

    if gen_u2: p = 2
    else: p = 1

    N = signal.shape[0] # signal length
    L = ul.shape[0] # number of basis functions
    K = int((N - 1 - p)/2.)


    assert (2*K + p < N), "choose higher p, need more signal points"

    if gen_u2:
        assert p >= 2, "choose higher p, U2 requires p=2"
    
    ul_inv = 1./ul
    ul_invK = zeros_like(ul_inv, dtype=complex128)
    for l in xrange(L):
        ul_invK[l] = accurate_pow(ul_inv[l], K)
    
    ul_invk = ones_like(ul_inv)
    g0 = zeros((L,), dtype=complex128)
    g0_K = zeros_like(g0)
    
    D0 = zeros_like(g0)
    U0 = zeros((L, L), dtype=complex128)
    U1 = zeros_like(U0)

    if gen_u2:
        g1 = zeros((L,), dtype=complex128)
        g1_K = zeros_like(g0)
        U2 = zeros_like(U0)
    
    for k in xrange(K + 1): # iterate over signal halves and accumulate G0, G0_M
        for l in xrange(L):
            g0[l] += ul_invk[l]*signal[k]
            g0_K[l] += ul_invk[l]*signal[K + 1 + k]
            
            D0[l] += (k + 1)*signal[k]*ul_invk[l] \
                + (K - k)*signal[k + K + 1]*ul_invk[l]*ul_inv[l]*ul_invK[l]

            if gen_u2:
                g1[l] += ul_invk[l]*signal[k+1]
                g1_K[l] += ul_invk[l]*signal[K + 1 + k + 1]

            if k % NPOW == NPOW-1:
                ul_invk[l] = accurate_pow(ul_inv[l], k + 1)
            else:
                ul_invk[l] *= ul_inv[l]
    
    for l in xrange(L):
        for lp in xrange(l):
            U0[l, lp] = 1./(ul[l] - ul[lp]) * (ul[l]*g0[lp] - ul[lp]*g0[l] 
                + ul_invK[lp]*g0_K[l] - ul_invK[l]*g0_K[lp])
            U0[lp, l] = U0[l, lp]
            
            U1[l, lp] = 1./(ul[l] - ul[lp]) * (ul[l]*ul[lp]*g0[lp] - ul[lp]*ul[l]*g0[l] 
                - ul_invK[l]*ul[lp]*g0_K[lp] + ul_invK[lp]*ul[l]*g0_K[l])
            U1[lp, l] = U1[l, lp]

            if gen_u2:
                U2[l, lp] = 0.5*((ul[l] + ul[lp])*U1[l, lp] - ul[l]*g1[lp] \
                    - ul[lp]*g1[l] + ul_invK[l]*g1_K[lp] + ul_invK[lp]*g1_K[l])
                U2[lp, l] = U2[l, lp]
            
        U0[l, l] = D0[l]
        U1[l, l] = D0[l]*ul[l] - ul[l]*g0[l] + ul_invK[l]*g0_K[l]

        if gen_u2:
            U2[l, l] = ul[l]*U1[l, l] - ul[l]*g1[l] + ul_invK[l]*g1_K[l]
        
    if gen_u2:
        return (U0, U1, g0, U2)
    else:
        return (U0, U1, g0, g0_K, D0)

def gen_amplitudes(B, g0):
    assert B.shape[0] == g0.shape[0], 'incorrect shape for B'

    L = g0.shape[0]
    N = B.shape[1]
    
    amps = zeros((N,), dtype=complex128)
    
    for n in xrange(N):
        amps[n] = dot(B[:, n], g0)**2
    return amps

def generate_rel_err(eigs, vr, u2):
    assert eigs.shape[0] == vr.shape[0], \
        'number of eigvals and eigvecs must match'
    assert vr.shape[0] == u2.shape[0], \
        'number of eigvecs and U2 rows must be compatible'

    rel_err = zeros((eigs.shape[0],), dtype=float64)

    for l in xrange(eigs.shape[0]):
        VU2V = dot(vr[:, l], dot(u2, vr[:, l]))
        rel_err[l] = abs(0.5*log(VU2V / eigs[l]**2))/abs(eigs[l])

    return rel_err

def make_basis(fmin, fmax, L=None, T=None, dt=1.0, tweak=0.1):
    '''generate a basis list with equally-spaced eigenfrequencies in the
    interval [fmin, fmax]

    tweak parameter notes the fractional increase in the number of basis
    functions from the standard estimate $(N dt / 4 \pi) * (fmax-fmin)$
    '''
    assert fmax > fmin, "fmax must be larger than fmin"
    assert (2*pi*fmax*dt < pi), "(2*Nyquist) criterion violated" # not quite true

    rho = 1.0 + tweak # spectral density. 1 corresponds to Nyquist

    if L is None:
        # calculate average density from time-bandwidth product
        assert T is not None, "must specify total length for L calculation"

        N = int(T/dt)
        frac_freq = (0.5/T*2*pi*(fmax-fmin)*dt)
        L = int(rho*N*frac_freq)
        '''
        print("N = ", N, " sample points")
        print("Using ", L, " basis functions")
        print("That's ", L/(frac_freq*N), " spectral density")
        '''
        assert L > 1, "change frequency window to increase number of basis \
        functions"

    ul = exp(-1j*2*pi*dt*linspace(fmin, fmax, L))
    return ul

def solve_fdm(ul, signal, amplitude_threshold=1e-8, relerr_threshold=1e-8, 
        gen_u2=False):

    # TODO: actually these numbers are related to condition number 
    # from ZGGEV, so I can check spuriousness that way too
    SMALL = 1e-14
    LARGE = 1e+14

    if gen_u2:
        U0, U1, g0, U2 = generate_u(ul, signal, gen_u2=True)
    else:
        U0, U1, g0 = generate_u(ul, signal, gen_u2=False)

    alpha,beta,vl,vr,work,info = zggev(U1, U0, compute_vl=False, compute_vr=True)

    if info != 0:
        print("zggev failed with code ", info)
    
    # remove zero or infinite eigenvalues
    idx = where((abs(alpha) < LARGE) & (abs(beta) > SMALL))[0]
    eigs = (log(alpha[idx]) - log(beta[idx]))/(-1j*2*pi)

    if idx.shape[0] != alpha.shape[0]:
        print("there were {:d} singular values removed" \
            .format(alpha.shape[0] - idx.shape[0]))
    
    vr = vr[:, idx]

    for l in xrange(vr.shape[1]):
        # rescale vectors
        norm = dot(vr[:, l], dot(U0, vr[:, l]))
        vr[:, l] *= 1.0/sqrt(norm)

    #assert abs(dot(dot(vr[:, 0], U0), vr[:, 0]) - 1) < 1e-3, "incorrect normalization"
    if gen_u2:
        errs = generate_rel_err(alpha[idx]/beta[idx], vr, U2) 

    A = gen_amplitudes(vr, g0)

    if gen_u2:
        idx = where((abs(A) > amplitude_threshold) & (errs < relerr_threshold))
        return eigs[idx].squeeze(), A[idx].squeeze(), errs[idx]
    else:
        idx = where((abs(A) > amplitude_threshold))
        return eigs[idx].squeeze(), A[idx].squeeze()

def make_lorenzian(N, amplitudes, poles, dt=1.):
    n = arange(N)
    signal = zeros((N,), dtype=complex128)

    try:
        for a,f in it.izip(amplitudes, poles):
            signal += a*exp(-1j*2*pi*f*dt*n)
    except TypeError:
        signal = amplitudes*exp(-1j*2*pi*poles*dt*n)
    
    return signal

