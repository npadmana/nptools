#ifndef SPDATA_H
#define SPDATA_H

#include <Eigen/Dense>

class SpData
{
public:
    double z;
    Eigen::VectorXd wave, flux, ivar;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_NOTHROW; ///??????

    // Read in an spData object from a file
    SpData(std::string fn);

    writeto(std::string fn);
    template <class F> writeto(std::string fn, F func);

    // Compute chi2
    template <class F> double chi2(F func) {;} //?????
};

#endif // SPDATA_H
