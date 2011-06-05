#include <gtest/gtest.h>
#include "linelist.h"
#include <stdexcept>


TEST(LineList, WaveCount) {
    LineList ll;
    EXPECT_EQ(ll.nlines("NoSuchLine"), 0);
    EXPECT_THROW(ll.wave("NoSuchLine"), std::runtime_error);

    EXPECT_EQ(ll.nlines("CIV"), 1);
    EXPECT_EQ(ll.wave("CIV").size(), 1);
    EXPECT_EQ(ll.wave("CIV")[0], 1545.86);

    EXPECT_EQ(ll.nlines("CIII"), 1);
    EXPECT_EQ(ll.wave("CIII").size(), 1);
    EXPECT_EQ(ll.wave("CIII")[0], 1908.27);

    EXPECT_EQ(ll.nlines("MgII"), 1);
    EXPECT_EQ(ll.wave("MgII").size(), 1);
    EXPECT_EQ(ll.wave("MgII")[0], 2800.32);



}
