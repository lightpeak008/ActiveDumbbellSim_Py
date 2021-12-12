
import cython
import numpy as np
cimport numpy as np
np.import_array()

from libc.math cimport round, sqrt, ceil, floor

DOUBLE = np.float64
INTEGER = np.int32


cdef class LJBox:
    cdef double epsilon, sigma
    cdef double r_cut, r_cut2
    cdef double e_cut
    cdef int n_part, n_dim
    cdef double[:] box_shape

    cdef double[:] rij
    cdef double[:] cell_shape
    cdef int[:] cell_num
    cdef int[:] cell_tot

    cdef int[:] cell_heads
    cdef int[:] cell_links
    cdef int[:, :] iarr2lat
    cdef int[:] iarr_part
    cdef int[:] iarr_dim

    def __init__(self, double epsion, double sigma, double r_cut,
                 int n_dim, double[:] box_shape, int n_part):
        self.epsilon, self.sigma = epsion, sigma
        self.r_cut = r_cut
        self.r_cut2 = self.r_cut * self.r_cut

        self.n_dim, self.n_part = n_dim, n_part
        self.box_shape = box_shape

        self.rij = np.zeros(self.n_dim, dtype=DOUBLE)

        cdef int i_dim
        self.cell_num = np.zeros(n_dim, dtype=INTEGER)
        self.cell_shape = np.zeros(n_dim, dtype=DOUBLE)
        self.cell_tot = 1
        for i_dim in range(self.n_dim):
            self.cell_num[i_dim] = ceil(self.box_shape[i_dim] / self.r_cut)
            self.cell_tot += self.cell_num[i_dim]
            self.cell_shape[i_dim] = self.box_shape[i_dim] / self.cell_num[i_dim]

        self.cell_heads = np.zeros(self.n_part, dtype=INTEGER)
        self.cell_links = np.zeros(self.n_part, dtype=INTEGER)

        self.iarr2lat = np.zeros((self.cell_tot, self.n_dim), dtype=INTEGER)
        cdef int iarr, fac
        fac = self.cell_tot
        for iarr in range(self.cell_tot):
            for i_dim in range(self.n_dim - 1, -1, -1):
                fac /= self.cell_num[i_dim]
                self.iarr2lat[iarr, i_dim] = iarr // fac
                iarr = iarr % fac

        self.iarr_part = np.zeros(self.n_part, dtype=INTEGER)
        self.iarr_dim = np.zeros(self.n_dim, dtype=INTEGER)

        cdef double r2, r6
        r2 = self.r_cut2 / (self.sigma * self.sigma)
        r6 = r2 * r2 * r2
        self.e_cut = 4 * self.epsilon * (1 / (r6 * r6) - 1 / r6)

    cdef int lat2arr(self, int[:] lat):
        cdef int i_dim, fac, iarr
        fac = 1
        iarr = 0
        for i_dim in range(self.n_dim):
            iarr += fac * lat[i_dim]
            fac *= self.cell_num[i_dim]
        return iarr

    cdef void gen_cell(self, double[:, ::1] pos):
        cdef int[:] iarr_sort

        cdef int i_part, i_dim
        for i_part in range(self.n_part):
            for i_dim in range(self.n_dim):
                self.iarr_dim[i_dim] = floor(pos[i_part, i_dim] / self.cell_shape[i_dim])
            self.iarr_part[i_part] = self.lat2arr(self.iarr_dim)

        iarr_sort = np.argsort(self.iarr_part)

        cdef int i_cell = 0
        cdef int start = 0
        cdef int end = 0
        cdef int j_part
        for i_part in range(self.n_part):
            if i_cell == self.iarr_part[iarr_sort[i_part]]:
                end = i_part
            else:
                self.cell_heads[i_cell] = iarr_sort[end]
                for j_part in range(end, start):
                    self.cell_links[iarr_sort[j_part]] = self.cell_links[iarr_sort[j_part-1]]
                self.cell_links[iarr_sort[start]] = -1

                start = i_part
                end = i_part
                i_cell = self.iarr_part[iarr_sort[i_part]]

    cdef