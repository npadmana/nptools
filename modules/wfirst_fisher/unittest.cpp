#include <gtest/gtest.h>
#include <stdexcept>
#include <cmath>
#include "zspace.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Test case based on the example in the code.
TEST(Zspace, zspace1) {
    VectorXd nbar(1), bias(1);
    nbar << 5.e-4; bias << 1.5;

    MatrixXd cov(2,2);

    // Execute the function
    double f = pow(0.25, 0.6);
    zspace_mbias_pk(nbar, 0.8, bias, f, 0.0, 1.e9, 0.1, cov);

    EXPECT_NEAR(0.0111924, sqrt(cov(0,0)), 1.e-5);
    EXPECT_NEAR(0.0272365, sqrt(cov(1,1)), 1.e-5);
}
