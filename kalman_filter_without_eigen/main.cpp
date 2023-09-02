#include <iostream>
#include <vector>
#include <tuple>
#include "matrix"

// Kalman filter function
std::pair<Matrix, Matrix> kalman_filter(Matrix& x,  Matrix& P,  Matrix& u,
                                        Matrix& F,  Matrix& H,  Matrix& R,  Matrix& I,
                                       const std::vector<double>& measurements) {
    for (size_t n = 0; n < measurements.size(); n++) {
        // Measurement Update
        Matrix Z(1, 1);
        Z(0, 0) = measurements[n];

        Matrix y = Z - H*x; // Measurement residual
        Matrix S = (H * P * H.T()) + R; // Measurement covariance

        // Calculate Kalman gain using MatrixInverse function
        Matrix K = P * H.T() * S.I(); // Kalman gain

        x = x + K * y; // Update state estimate
        P = P - K * (H * P); // Update covariance matrix

        // Prediction
        x = F * x + u;      // Predict state
        P = F * P * F.T();    // Predict covariance
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
