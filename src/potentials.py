import numpy as np
from scipy.spatial.distance import pdist, squareform


class LJPotential:
    r_small = 1e-6

    def __init__(self, epsilon, sigma, r_cut, r_skin,
                 box_shape, n_part, n_dim):
        self.eps = epsilon
        self.sig = sigma

        self.r_cut, self.r_skin = r_cut, r_skin

        self.box_shape = box_shape
        self.n_part, self.n_dim = n_part, n_dim

        self.verlet_list = None

    def gen_verlet_list(self, pos):


    def compute_int(self, pos):
        rij = pos.reshape(1, self.n_part, self.n_dim) - pos.reshape(self.n_part, 1, self.n_dim)

