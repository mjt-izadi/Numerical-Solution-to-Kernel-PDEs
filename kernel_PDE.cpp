#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

class KernelParabolic {
public:
    // Backstepping kernel for parabolic PDE
    KernelParabolic(double lam, double c, int N) : lam(lam), c(c), N(N) {
        h = 1.0 / (N - 1);
        xi = ArrayXd::LinSpaced(2 * N - 1, 0, 2);
        eta = ArrayXd::LinSpaced(N, 0, 1);
        generateIndices();
    }

    void calculateKernelXiEta() {
        const double const_val = (lam + c) / 4;
        A = MatrixXd::Zero(N * N, N * N);
        b = VectorXd::Zero(N * N);

        for (int j = 0; j < N; ++j) {
            for (int i = j; i < 2 * N - j - 1; ++i) {
                int id = indices(i, j);
                A(id, id) -= 1.0;
                b(id) -= const_val * (-xi(i) + eta(j));

                for (int m = j; m < i; ++m) {
                    for (int n = 0; n < j; ++n) {
                        A(id, indices(m, n)) += const_val * h * h / 4;
                        A(id, indices(m, n + 1)) += const_val * h * h / 4;
                        A(id, indices(m + 1, n)) += const_val * h * h / 4;
                        A(id, indices(m + 1, n + 1)) += const_val * h * h / 4;
                    }
                }
            }
        }
        kernel = A.fullPivLu().solve(b);
    }

    void calculateGain() {
        y = VectorXd::Zero(N);
        gain = VectorXd::Zero(N);
        for (int n = 0; n < N; ++n) {
            y(n) = (xi(N + n - 1) - eta(N - n - 1)) / 2;
            int idx = indices(N + n - 1, N - n - 1);
            gain(n) = kernel(idx);
        }
    }

    void printResults() {
        std::cout << "y: " << std::endl << y.transpose() << std::endl;
        std::cout << "gain: " << std::endl << gain.transpose() << std::endl;
    }

private:
    double lam;
    double c;
    int N;
    double h;
    ArrayXd xi;
    ArrayXd eta;
    MatrixXi indices;
    MatrixXd A;
    VectorXd b;
    VectorXd kernel;
    VectorXd y;
    VectorXd gain;

    void generateIndices() {
        indices = MatrixXi::Zero(2 * N - 1, N);
        for (int j = 0; j < N; ++j) {
            for (int i = j; i < 2 * N - j - 1; ++i) {
                indices(i, j) = j * (2 * N - j - 1) + i;
            }
        }
    }
};

int main() {
    KernelParabolic k1(10.0, 0.0, 25);
    k1.calculateKernelXiEta();
    k1.calculateGain();
    k1.printResults();
    return 0;
}
