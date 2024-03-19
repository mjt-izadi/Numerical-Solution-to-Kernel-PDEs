import numpy as np
import matplotlib.pyplot as plt

class Kernel_parabolic:
    # Backstepping kernel for parabolic PDE
    def __init__(self, lam, c, N):
        self.lam = lam
        self.c = c
        self.N = N
        self.h = 1 / (self.N - 1)
        self.xi = np.linspace(0, 2, 2 * N - 1)
        self.eta = np.linspace(0, 1, N)
        self.generate_indices()

    def generate_indices(self):
        self.indices = np.zeros((2 * self.N - 1, self.N), dtype=int)
        for j in range(self.N):
            for i in range(j, 2 * self.N - j - 1):
                self.indices[i, j] = j * (2 * self.N - j - 1) + i

    def calculate_kernel_xieta(self):
        const = (self.lam + self.c) / 4
        A = np.zeros((self.N ** 2, self.N ** 2))
        b = np.zeros(self.N ** 2)
        for j in range(self.N):
            for i in range(j, 2 * self.N - j - 1):
                A[self.indices[i, j], self.indices[i, j]] += -1.0
                b[self.indices[i, j]] += -const * (-self.xi[i] + self.eta[j])
                for m in range(j, i):
                    for n in range(j):
                        A[self.indices[i, j], self.indices[m, n]] += const * self.h ** 2 / 4
                        A[self.indices[i, j], self.indices[m, n + 1]] += const * self.h ** 2 / 4
                        A[self.indices[i, j], self.indices[m + 1, n]] += const * self.h ** 2 / 4
                        A[self.indices[i, j], self.indices[m + 1, n + 1]] += const * self.h ** 2 / 4
        self.kernel = np.linalg.solve(A, b)

    def calculate_gain(self):
        ij = np.array([[self.N + n - 1, self.N - n - 1] for n in range(self.N)])
        self.y = (self.xi[ij[:, 0]] - self.eta[ij[:, 1]]) / 2
        gain_indices = [self.indices[i, j] for i, j in ij]
        self.gain = self.kernel[gain_indices]            


if __name__ == "__main__":
    N = 50
    k1 = Kernel_parabolic(25.0, 0.0, N)
    k1.calculate_kernel_xieta()
    k1.calculate_gain()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(k1.y, k1.gain)
    plt.show()