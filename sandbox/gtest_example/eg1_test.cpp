#include "eg1.h"
#include "gtest/gtest.h"

TEST(Simple, TimesTwo) {
    EXPECT_EQ(timestwo(2), 4);
}
