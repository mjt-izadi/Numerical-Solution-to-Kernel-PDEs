import numpy as np
import matplotlib.pyplot as plt

class Kernel_parabolic:
    # Backstepping kernel for parabolic PDE
    lam = 10.0
    c = 0.0

    def __init__(self, N):
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
                        A[self.indices[i, j], self.indices[m, n]] += const * self.h / 4
                        A[self.indices[i, j], self.indices[m, n + 1]] += const * self.h / 4
                        A[self.indices[i, j], self.indices[m + 1, n]] += const * self.h / 4
                        A[self.indices[i, j], self.indices[m + 1, n + 1]] += const * self.h / 4
        self.kernel = np.linalg.solve(A, b)

if __name__ == "__main__":
    N = 50
    k1 = Kernel_parabolic(N)
    k1.calculate_kernel_xieta()

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(k1.kernel)
    plt.show()