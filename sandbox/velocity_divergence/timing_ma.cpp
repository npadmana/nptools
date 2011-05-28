//#define BOOST_DISABLE_ASSERTS
#include <ctime>
#include <cstdio>
#include <boost/lambda/lambda.hpp>
#include "multi_array_wrap.h"

using namespace ma;

int main(int argc, char* argv[])
{
    const int X_SIZE = 2000;
    const int Y_SIZE = 2000;
    const int ITERATIONS = 50;
    unsigned int startTime = 0;
    unsigned int endTime = 0;

    // Create the boost array
    typedef boost::multi_array<double, 2> ImageArrayType;
    ImageArrayType boostMatrix(boost::extents[X_SIZE][Y_SIZE]);
    MA<ImageArrayType> test(boost::extents[X_SIZE][Y_SIZE]);

    // Create the native array
    double *nativeMatrix = new double [X_SIZE * Y_SIZE];

    //------------------Measure MA----------------------------------------------
    startTime = clock();
    for (int i = 0; i < ITERATIONS; ++i)
    {
        multi_for(test, _1=2.345);
    }
    endTime = clock();
    printf("[MA] Elapsed time: %6.3f seconds\n", (endTime - startTime) / (double)CLOCKS_PER_SEC);

    //------------------Measure boost----------------------------------------------
    startTime = clock();
    for (int i = 0; i < ITERATIONS; ++i)
    {
        for (int y = 0; y < Y_SIZE; ++y)
        {
            for (int x = 0; x < X_SIZE; ++x)
            {
                boostMatrix[y][x] = 2.345;
            }
        }
    }
    endTime = clock();
    printf("[Boost] Elapsed time: %6.3f seconds\n", (endTime - startTime) / (double)CLOCKS_PER_SEC);

    //------------------Measure native-----------------------------------------------
    startTime = clock();
    for (int i = 0; i < ITERATIONS; ++i)
    {
        for (int y = 0; y < Y_SIZE; ++y)
        {
            for (int x = 0; x < X_SIZE; ++x)
            {
                nativeMatrix[x + (y * X_SIZE)] = 2.345;
            }
        }
    }
    endTime = clock();
    printf("[Native]Elapsed time: %6.3f seconds\n", (endTime - startTime) / (double)CLOCKS_PER_SEC);

    return 0;
}

