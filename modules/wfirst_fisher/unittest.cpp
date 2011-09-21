#include <gtest/gtest.h>
#include <stdexcept>
#include <cmath>
#include "zspace.h"
#include "bao_forecast.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Put in some basic tests of cosmological terms
TEST(Cosmo, comdis) {
    detf fidcosmo = detf_fiducial();
    EXPECT_NEAR(0.44785/0.719, comdis(z2a(0.5), fidcosmo), 1.e-3); // comdis returns distances without h^-1 in them
    EXPECT_NEAR(0.63208/0.719, comdis(z2a(0.75), fidcosmo), 1.e-3); // comdis returns distances without h^-1 in them
    EXPECT_NEAR(0.79247/0.719, comdis(z2a(1.0), fidcosmo), 1.e-3); // comdis returns distances without h^-1 in them
}

TEST(Cosmo, volume) {
    detf fidcosmo = detf_fiducial();
    EXPECT_NEAR(11.18, shellVol_Gpc_h(z2a(0.5), z2a(1.0), 10000.0, fidcosmo), 1.e-2);
    EXPECT_NEAR(11.18/2.0, shellVol_Gpc_h(z2a(0.5), z2a(1.0), 5000.0, fidcosmo), 1.e-2); // Scaling test
}

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

// Test sigma_perp against some old values I had
TEST(BAOSigma, sperp_z) {
    EXPECT_NEAR(5.766, Sigma_perp(z2a(0.75)), 1.e-1);
    EXPECT_NEAR(4.5, Sigma_perp(z2a(1.35)), 1.e-1);
    EXPECT_NEAR(4.04, Sigma_perp(z2a(1.65)), 1.e-1);
}

// Test sigma_par against some old values I had
// They were calculated somewhat differently (Om^0.6, instead of Om^0.55)
// so that explains the low accuracy.
TEST(BAOSigma, spar_z) {
    EXPECT_NEAR(10.128, Sigma_par(z2a(0.75)), 3.e-1);
    EXPECT_NEAR(8.456, Sigma_par(z2a(1.35)), 3.e-1);
    EXPECT_NEAR(7.718, Sigma_par(z2a(1.65)), 3.e-1);
}


// Now test the BAO forecast code
// The original code that produced these results is in originals
TEST(BAOFisher, test1) {
    Matrix2d tmp = bao_forecast(1.e-4, 0.8, 5.0, 10.0, 0.0, 0.0, 1.0);
    EXPECT_NEAR(5.729, sqrt(tmp(0,0)), 1.e-3);
    EXPECT_NEAR(12.467, sqrt(tmp(1,1)), 1.e-3);
    EXPECT_NEAR(0.40529, tmp(0,1)/sqrt(tmp(0,0)*tmp(1,1)), 1.e-3);
    EXPECT_NEAR(0.40529, tmp(1,0)/sqrt(tmp(0,0)*tmp(1,1)), 1.e-3);
}

TEST(BAOFisher, test2) {
    Matrix2d tmp = bao_forecast(1.e-4, 0.8, 5.0, 10.0, 0.0, 1.0, 1.0);
    EXPECT_NEAR(5.0316, sqrt(tmp(0,0)), 1.e-3);
    EXPECT_NEAR(7.5236, sqrt(tmp(1,1)), 1.e-3);
    EXPECT_NEAR(0.4317, tmp(0,1)/sqrt(tmp(0,0)*tmp(1,1)), 1.e-3);
    EXPECT_NEAR(0.4317, tmp(1,0)/sqrt(tmp(0,0)*tmp(1,1)), 1.e-3);
}


// A case from my old python code
TEST(BAOFisher, test3) {
    Matrix2d tmp = bao_forecast_shell(10.95e-4, 0.608, 0.0, 0.681, 1.0 ,1.1, 24000.0, 0.5, false);
    EXPECT_NEAR(0.548, sqrt(tmp(0,0)), 3.e-2);
    EXPECT_NEAR(0.861, sqrt(tmp(1,1)), 3.e-2);
    EXPECT_NEAR(0.419, tmp(0,1)/sqrt(tmp(0,0)*tmp(1,1)), 1.e-3);
    EXPECT_NEAR(0.419, tmp(1,0)/sqrt(tmp(0,0)*tmp(1,1)), 1.e-3);
}

