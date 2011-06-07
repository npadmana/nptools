#include <gtest/gtest.h>
#include <stdexcept>
#include "forecast.h"

using namespace std;

TEST(AZ, a2z) {
    EXPECT_NEAR(0.0, a2z(1.0), 1.e-7);
    EXPECT_NEAR(1.0, a2z(0.5), 1.e-7);
}

TEST(AZ, z2a) {
    EXPECT_NEAR(1.0, z2a(0.0), 1.e-7);
    EXPECT_NEAR(0.5, z2a(1.0), 1.e-7);
}


TEST(FoMSWGTest, RedshiftZero) {
    Fswg cc = Fiducial();
    EXPECT_NEAR(0.719, hubble(1.0, cc), 1.e-4);
    EXPECT_NEAR(0.2565, Omega_M0(cc), 1.e-4);
    EXPECT_NEAR(0.2565, Omega_M(1.0, cc), 1.e-4);
    EXPECT_NEAR(Omega_M0(cc), Omega_M(1.0, cc), 1.e-4);
}


TEST(FoMSWGTest, HighZ) {
    Fswg cc = Fiducial();
    EXPECT_NEAR(1.0, Omega_M(0.001, cc), 1.e-4);
}

TEST(FoMSWGTest, Growth) {
    Fswg cc = Fiducial();
    EXPECT_NEAR(1.0, fgrowth(0.001, cc), 1.e-4);
    // Growth factor at high redshift is ~ a
    EXPECT_NEAR(0.001, growth(0.001, cc), 1.e-4);
}

TEST(Forecast, SigmaPP) {
    Fswg cc = Fiducial();
    const double s0 = 11.022222;
    EXPECT_NEAR(s0*growth(1.0, cc), SigmaPerp(1.0), 1.e-3);
    EXPECT_NEAR(s0*(1+fgrowth(1.0, cc))*growth(1.0, cc), SigmaPar(1.0), 1.e-3);
}
