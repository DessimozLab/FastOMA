"""
OMA HOG distance fitting
Fits distances from a subset of trees to OMA HOGs.

(C) 2025 Comparative Genomics Group <contact@omabrowser.org>

This file is part of hog_distances.

hog_distances is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

hog_distances is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with hog_distances.  If not, see <http://www.gnu.org/licenses/>.
"""

from numba import jit, prange
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsmr
import numpy as np


def nngs(A, d): #, edges_fixed, fit_missing_free):  # converge_type=2, initial=1):
    #func = _nngs1 if converge_type == 1 else _nngs2
    func = _nngs2

    H = (A.T@A)
    if type(H) is not csr_matrix:
        H = H.tocsr()
    D = H.diagonal()
    ii = H.indptr
    jj = H.indices
    xx = H.data

    f = -1.0*(A.T@d)

    #x0 = np.zeros_like(f)
    return func(ii, jj, xx, D, f, H.shape[0])#, x0)

    #mask = np.zeros_like(x0, dtype=np.bool_)
    #for (i, e_i) in enumerate(edges_fixed):
    #    x0[e_i] = d[i]
    #    if not fit_missing_free:
    #        mask[e_i] = True



    #mask = edges_fixed

    #return func(ii, jj, xx, D, f, H.shape[0], x0, mask)

    #if initial == 1:
    #    # first try with unconstrained soln
    #    x0 = lsmr(A, d)[0]  # starting guess - unconstrained least squares
    #    ret = func(ii, jj, xx, D, f, H.shape[0], x0)
    #
    #    r = ret[2]
    #    if r < 1e-6:
    #        return ret
    #    else:
    #        x0s = np.zeros_like(x0)
    #        ret1 = func(ii, jj, xx, D, f, H.shape[0], x0s)
    #
    #        r1 = ret1[2]
    #        return (ret1 if r1 < r else ret)
    #elif initial == 2:
    #    x0 = lsmr(A, d)[0]
    #    return func(ii, jj, xx, D, f, H.shape[0], x0)
    #elif initial == 3:
    #    x0 = np.zeros_like(f)
    #    return func(ii, jj, xx, D, f, H.shape[0], x0)


#@jit(nopython=True, fastmath=True, nogil=True, parallel=True)
#def _nngs1(ii, jj, xx, D, f, n, x):
#    EPS = 1e-6
#
#    mu = np.copy(f)
#    #x = np.zeros_like(f)
#    #x1 = np.zeros_like(f)
#    x1 = np.zeros_like(x)
#
#    it = 0
#    converge = False
#    while not converge:
#        it += 1
#        for i in range(n):
#            x1[i] = max(0.0, (x[i] - (mu[i] / D[i])))
#            s = x1[i] - x[i]
#            for j in prange(ii[i], ii[i+1]):
#                mu[jj[j]] += (s * xx[j])
#
#        converge = test_convergence(x, x1, EPS)
#
#        # Swap array
#        tx = x
#        x = x1
#        x1 = tx
#
#    m = np.inf
#    z = 0
#    for i in range(n):
#        t = 0
#        for j in range(ii[i], ii[i+1]):
#            t += xx[j] * x[jj[j]]
#        z += (t * x[i])
#        t += f[i]
#        if t < m:
#            m = t
#    r = (z + np.dot(x, f) - (np.sum(x)*m))
#    return (x, it, r)


#@jit(nopython=True, fastmath=True, nogil=True, parallel=True)
#def test_convergence(x, y, eps):
#    # convergence test -- usual numpy does not work here in numba...
#    z = np.abs(x - y)
#    c = 0
#    for i in prange(z.shape[0]):
#        c += (1 if z[i] >= eps else 0)
#
#    return (c == 0)


@jit(nopython=True, fastmath=True, nogil=True, parallel=True)
def _nngs2(ii, jj, xx, D, f, n):#, x):#, mask):
    #DELTA = 1e-9
    DELTA = 1e-3
    ITER_CHECK = min(n // 2, 500)
    MAXITERS = int(1e3)

    mu = np.copy(f)
    x = np.zeros_like(f)

    it = 0
    converge = False
    r1 = np.inf
    while not converge and it < MAXITERS:
        it += 1
        for i in range(n):
            t = max(0, x[i] - (mu[i] / D[i]))
            s = t - x[i]
            #if not mask[i]:
            x[i] = t
            for j in prange(ii[i], ii[i+1]):
                mu[jj[j]] += (s * xx[j])

        if it % ITER_CHECK > 0:
            pass
        else:
            # Convergence test
            m = np.inf
            z = 0
            for i in range(n):
                t = 0
                for j in prange(ii[i], ii[i+1]):
                    t += xx[j] * x[jj[j]]
                z += (t * x[i])
                t += f[i]
                if t < m:
                    m = t
            r = (z + np.dot(x, f) - (np.sum(x)*m))
            if r < DELTA or r == r1:
                converge = True
            r1 = r

    return (x, it, r)
