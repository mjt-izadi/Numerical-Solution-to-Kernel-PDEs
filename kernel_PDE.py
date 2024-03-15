import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle

class kernel_TVD:
    tau = 1.0
    alpha = 4.0

    h = 0.001
    n_t = 10000

    def __init__(self, N, set_number):
        self.N = N
        self.set_number = set_number
        # self.xi1 = np.arange(0, 2 * N - 1)
        # self.eta1 = np.arange(0, N)
        self.beta = np.sqrt(self.alpha / self.tau)
        self.define_length()
        self.calculate_functions()
        self.generate_indices()

    def define_length(self):
        self.t = np.linspace(0, self.n_t * self.h, self.n_t)
        ind = np.nonzero(self.t <= 8.5)
        c = np.array([-0.00551028, -0.00068658])
        self.l = np.zeros_like(self.t)
        self.l_dot = np.zeros_like(self.t)
        self.l_ddot = np.zeros_like(self.t)
        self.l[ind] = c[0] * (self.t[ind] - 8.5) ** 3 +\
             c[1] * (self.t[ind] - 8.5) ** 4
        self.l += 1.2
        self.l_dot[ind] = 3 * c[0] * (self.t[ind] - 8.5) ** 2 +\
             4 * c[1] * (self.t[ind] - 8.5) ** 3
        self.l_ddot[ind] = 6 * c[0] * (self.t[ind] - 8.5) +\
             12 * c[1] * (self.t[ind] - 8.5) ** 2
        # self.l = np.sin(np.pi / 13 * self.t) + 1
        # self.l_dot =  np.pi / 13 * np.cos(np.pi / 13 * self.t)
        # self.l_ddot =  -np.pi * np.pi / 13 / 13 * np.sin(np.pi / 13 * self.t)

    def generate_indices(self):
        self.k1_ind = np.zeros((2 * self.N - 1, self.N), dtype=int)
        self.k2_ind = np.zeros_like(self.k1_ind)
        for j in range(self.N - 2 + 1):
            for i in range(j + 1, 2 * self.N - j - 2 + 1):
                self.k1_ind[i, j] = j * (2 * self.N - j - 2) + i - 1
        for j in range(1, self.N - 1 + 1):
            for i in range(j, 2 * self.N - j - 2 + 1):
                self.k2_ind[i, j] = self.N * (self.N - 1) +\
                     (j - 1) * (2 * self.N - j - 2) + i - 1
        for j in range(self.N - 1 + 1):
            self.k1_ind[j, j] = self.N * (2 * self.N - 3) + j + 1
        for i in range(2 * self.N - 2 + 1):
            self.k2_ind[i, 0] = 2 * self.N *(self.N - 1) + i + 1

    def calculate_functions(self):
        self.fcn = {"lambda1":
                    1 / 2 * (-self.l_dot / np.sqrt(self.alpha * self.tau) +\
                             1 / self.tau) *\
                             np.exp(self.l / np.sqrt(self.alpha * self.tau)),
                    "lambda1_dot":
                    1 / 2 * (-self.l_ddot / np.sqrt(self.alpha * self.tau) +\
                             self.l_dot / np.sqrt(self.alpha * self.tau) *\
                             (-self.l_dot / np.sqrt(self.alpha * self.tau) +\
                              1 / self.tau)) *\
                              np.exp(self.l / np.sqrt(self.alpha * self.tau)),
                    "lambda2":
                    1 / 2 * (self.l_dot / np.sqrt(self.alpha * self.tau) +\
                             1 / self.tau) *\
                             np.exp(-self.l / np.sqrt(self.alpha * self.tau)),
                    "lambda2_dot":
                    1 / 2 * (self.l_ddot / np.sqrt(self.alpha * self.tau) +\
                             -self.l_dot / np.sqrt(self.alpha * self.tau) *\
                             (self.l_dot / np.sqrt(self.alpha * self.tau) +\
                              1 / self.tau)) *\
                              np.exp(-self.l / np.sqrt(self.alpha * self.tau)),
                    "f1": -np.exp(-self.l / np.sqrt(self.alpha * self.tau)),
                    "f1_dot":
                    self.l_dot / np.sqrt(self.alpha * self.tau) *\
                    np.exp(-self.l / np.sqrt(self.alpha * self.tau)),
                    "f2": -np.exp(self.l / np.sqrt(self.alpha * self.tau)),
                    "f2_dot":
                    -self.l_dot / np.sqrt(self.alpha * self.tau) *\
                    np.exp(self.l / np.sqrt(self.alpha * self.tau))}

    def calculate_ODE_matrices(self, i):
        Dx = self.l[i] / (self.N - 1)
        if self.set_number == 1:
            beta = self.beta
            lambda1 = self.fcn["lambda1"][i]
            lambda1_dot = self.fcn["lambda1_dot"][i]
            lambda2 = self.fcn["lambda2"][i]
            f = self.fcn["f1"][i]
            f_dot = self.fcn["f1_dot"][i]
        elif self.set_number ==2:
            beta = -self.beta
            lambda1 = self.fcn["lambda2"][i]
            lambda1_dot = self.fcn["lambda2_dot"][i]
            lambda2 = self.fcn["lambda1"][i]
            f = self.fcn["f2"][i]
            f_dot = self.fcn["f2_dot"][i]
        
        const = Dx / 4 / beta
        C = np.zeros((2 * self.N ** 2 - 3 * self.N + 1,
                      2 * self.N ** 2 - 3 * self.N + 1))
        D = np.zeros_like(C)
        F = np.zeros(2 * self.N ** 2 - 3 * self.N + 1)
        for I in range(1, 2 * self.N - 2 + 1):
            C[self.k1_ind[I, 0], self.k1_ind[1, 0]] += const # 1x
            D[self.k1_ind[I, 0], self.k1_ind[I, 0]] += -1.0 # 1x
            F[self.k1_ind[I, 0]] += -1 / 2 / beta * lambda1 *\
                (f + const * f_dot + 2 * const * I * lambda2) +\
                    -1 / 2 / beta * const * f * lambda1_dot # 1x
            for i in range(1, I - 1 + 1):
                C[self.k1_ind[I, 0], self.k1_ind[i, 0]] += const # 1x
                C[self.k1_ind[I, 0], self.k1_ind[i + 1, 0]] += const # 1x
        for J in range(1, self.N - 1 + 1):
            for I in range(J + 1, 2 * self.N - J - 2 + 1):
                C[self.k1_ind[I, J], self.k2_ind[J, J]] += const * f # 2x
                C[self.k1_ind[I, J], self.k1_ind[J + 1, J]] += const # 2x
                C[self.k2_ind[I, J], self.k2_ind[I, 1]] += const # 4x
                D[self.k1_ind[I, J], self.k1_ind[I, J]] += -1.0 # 2x
                D[self.k1_ind[I, J], self.k2_ind[J, J]] += f + const * f_dot # 2x
                D[self.k2_ind[I, J], self.k2_ind[I, J]] += -1.0 # 4x
                F[self.k2_ind[I, J]] += - 1 / 2 / beta *\
                    (lambda1 + const * lambda1_dot) # 4
                for i in range(J + 1, I - 1 + 1):
                    C[self.k1_ind[I, J], self.k1_ind[i, J]] += const # 2x
                    C[self.k1_ind[I, J], self.k1_ind[i + 1, J]] += const # 2x
                for j in range(1, J - 1 + 1):
                    C[self.k2_ind[I, J], self.k2_ind[I, j]] += const # 4x
                    C[self.k2_ind[I, J], self.k2_ind[I, j + 1]] += const # 4x
                for j in range(J - 1 + 1):
                    D[self.k2_ind[I, J], self.k1_ind[I, j]] += const * lambda1 # 4x
                    D[self.k2_ind[I, J], self.k1_ind[I, j + 1]] += const * lambda1 # 4x
                for i in range(J, I - 1 + 1):
                    D[self.k1_ind[I, J], self.k2_ind[i, J]] += const * lambda2 # 2x
                    D[self.k1_ind[I, J], self.k2_ind[i + 1, J]] += const * lambda2 # 2x
            C[self.k2_ind[J, J], self.k2_ind[J, 1]] += const # 3x
            D[self.k2_ind[J, J], self.k1_ind[J, J - 1]] += lambda1 * const # 3x
            D[self.k2_ind[J, J], self.k2_ind[J, J]] += -1.0 + lambda1 * const * f # 3x
            F[self.k2_ind[J, J]] += - 1 / 2 / beta *\
                 (lambda1 + const * lambda1_dot) # 3x
            for j in range(1, J - 1 + 1):
                C[self.k2_ind[J, J], self.k2_ind[J, j]] += const # 3x
                C[self.k2_ind[J, J], self.k2_ind[J, j + 1]] += const # 3x
            for j in range(J - 2 + 1):
                D[self.k2_ind[J, J], self.k1_ind[J, j]] += lambda1 * const # 3x
                D[self.k2_ind[J, J], self.k1_ind[J, j + 1]] += lambda1 * const # 3x
        return(C, D, F)

    def integrate_ODEs(self):
        self.k = np.zeros((2 * self.N ** 2 - 3 * self.N + 1, self.n_t))
        C, D, F = self.calculate_ODE_matrices(0)
        self.k[:, 0] = np.linalg.solve(D, -F)
        for i in tqdm(range(self.n_t - 1)):
            C, D, F = self.calculate_ODE_matrices(i)
            self.k[:, i +1] = np.linalg.solve(C + self.h * D,
                                        np.matmul(C, self.k[:, i]) - self.h * F)
            # self.k[:, i + 1] = np.linalg.solve(-12 * C + 5 * self.h * D,
            #                                    (np.matmul(-12 * C - 8 * self.h * D, self.k[:, i]) +\
            #                                    np.matmul(5 * self.h * D, self.k[:, i - 1]) +\
            #                                    -12 * self.h * F))

if __name__ == "__main__":
    N = 30
    # kernel1 = kernel_TVD(N, 1)
    # kernel1.integrate_ODEs()
    # with open('k1_30.pickle', 'wb') as file:
    #     pickle.dump(kernel1, file)
    kernel2 = kernel_TVD(N, 2)
    kernel2.integrate_ODEs()
    with open('k2_30.pickle', 'wb') as file:
        pickle.dump(kernel2, file)

    k = np.transpose(kernel2.k[4*N-8:4*N-5,:])
    tt = np.linspace(0, kernel2.n_t * kernel2.h, int(kernel2.n_t/10))
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(tt, k[::10])
    print(k[::10])
    plt.show()