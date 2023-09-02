#include <iostream>
#include <vector>
#include <tuple>

// Define a Matrix class to represent matrices and perform basic matrix operations
class Matrix {
public:
    Matrix(int rows, int cols) : rows(rows), cols(cols), data(rows * cols) {}

    double& operator()(int i, int j) {
        return data[i * cols + j];
    }

    double operator()(int i, int j) const {
        return data[i * cols + j];
    }

    int Rows() const {
        return rows;
    }

    int Cols() const {
        return cols;
    }

private:
    int rows;
    int cols;
    std::vector<double> data;
};

// Function to perform matrix multiplication
Matrix MatrixMultiply(const Matrix& a, const Matrix& b) {
    if (a.Cols() != b.Rows()) {
        throw std::invalid_argument("Matrix dimensions are incompatible for multiplication");
    }

    int m = a.Rows();
    int n = b.Cols();
    int p = a.Cols();
    Matrix result(m, n);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double sum = 0.0;
            for (int k = 0; k < p; ++k) {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }

    return result;
}

// Function to calculate the inverse of a 2x2 matrix
Matrix MatrixInverse(const Matrix& A) {
    double determinant = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

    if (determinant == 0.0) {
        throw std::runtime_error("Matrix is singular, cannot compute its inverse.");
    }

    Matrix result(2, 2);
    result(0, 0) = A(1, 1) / determinant;
    result(0, 1) = -A(0, 1) / determinant;
    result(1, 0) = -A(1, 0) / determinant;
    result(1, 1) = A(0, 0) / determinant;

    return result;
}

// Function to transpose a matrix
Matrix MatrixTranspose(const Matrix& a) {
    int rows = a.Rows();
    int cols = a.Cols();
    Matrix result(cols, rows);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(j, i) = a(i, j);
        }
    }

    return result;
}

// Function to add two matrices
Matrix MatrixAdd(const Matrix& a, const Matrix& b) {
    if (a.Rows() != b.Rows() || a.Cols() != b.Cols()) {
        throw std::invalid_argument("Matrix dimensions are incompatible for addition");
    }

    int rows = a.Rows();
    int cols = a.Cols();
    Matrix result(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result(i, j) = a(i, j) + b(i, j);
        }
    }

    return result;
}

// Kalman filter function
std::pair<Matrix, Matrix> kalman_filter(const Matrix& x, const Matrix& P, const Matrix& u,
                                       const Matrix& F, const Matrix& H, const Matrix& R, const Matrix& I,
                                       const std::vector<double>& measurements) {
    for (size_t n = 0; n < measurements.size(); n++) {
        // Measurement Update
        Matrix Z(1, 1);
        Z(0, 0) = measurements[n];

        Matrix y = Z - MatrixMultiply(H, x); // Measurement residual

        Matrix S = MatrixAdd(MatrixMultiply(H, MatrixMultiply(P, MatrixTranspose(H))), R); // Measurement covariance

        // Calculate Kalman gain using MatrixInverse function
        Matrix K = MatrixMultiply(P, MatrixMultiply(MatrixTranspose(H), MatrixInverse(S))); // Kalman gain

        x = x + MatrixMultiply(K, y); // Update state estimate
        P = P - MatrixMultiply(K, MatrixMultiply(H, P)); // Update covariance matrix

        // Prediction
        x = MatrixMultiply(F, x) + u; // Predict state
        P = MatrixMultiply(F, MatrixMultiply(P, MatrixTranspose(F))); // Predict covariance
    }

    return std::make_pair(x, P);
}

int main() {
    // Initialize matrices
    Matrix x(2, 1); // Initial state (location and velocity)
    x(0, 0) = 0.0;
    x(1, 0) = 0.0;

    Matrix P(2, 2); // Initial Uncertainty
    P(0, 0) = 100.0;
    P(0, 1) = 0.0;
    P(1, 0) = 0.0;
    P(1, 1) = 100.0;

    Matrix u(2, 1); // External Motion
    u(0, 0) = 0.0;
    u(1, 0) = 0.0;

    Matrix F(2, 2); // Next State Function
    F(0, 0) = 1.0;
    F(0, 1) = 1.0;
    F(1, 0) = 0.0;
    F(1, 1) = 1.0;

    Matrix H(1, 2); // Measurement Function
    H(0, 0) = 1.0;
    H(0, 1) = 0.0;

    Matrix R(1, 1); // Measurement Uncertainty
    R(0, 0) = 1.0;

    Matrix I(2, 2); // Identity Matrix
    I(0, 0) = 1.0;
    I(0, 1) = 0.0;
    I(1, 0) = 0.0;
    I(1, 1) = 1.0;

    // Actual measurements data
    std::vector<double> measurements = { 1.0, 2.0, 3.0 };

    // Call the Kalman filter function with measurements data
    std::pair<Matrix, Matrix> result = kalman_filter(x, P, u, F, H, R, I, measurements);

    // Print the results
    std::cout << "x=\n" << result.first(0, 0) << "\n" << result.first(1, 0) << "\n";
    std::cout << "P=\n" << result.second(0, 0) << " " << result.second(0, 1) << "\n"
              << result.second(1, 0) << " " << result.second(1, 1) << "\n";

    return 0;
}
