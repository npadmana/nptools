/** Uitilities to handle Eigen matrices. */

#ifndef EIGEN_UTILS_H
#define EIGEN_UTILS_H

#include <eigen3/Eigen/Core>
#include <string>
#include <fstream>
#include <boost/format.hpp>

/**
* Write out an Eigen derived class to a file fn, with a boost::format string fmt.
* This assumes a file format with ii, jj, M(ii, jj)
*/
template <typename DerivedBase>
void writeMatrix_withIndices(std::string fn, std::string fmt, const Eigen::DenseBase<DerivedBase>& mat) {
    double nrow, ncol;
    nrow = mat.rows();
    ncol = mat.cols();

    // Open file
    std::ofstream ff(fn.c_str());

    if (ff.is_open()) {
        for (int ii=0; ii < nrow; ++ii)
            for (int jj = 0; jj < ncol; ++jj)
                ff << boost::format(fmt) % ii % jj % mat(ii, jj);
    } else {
        throw "Unable to open file\n";
    }

    // Close file
    ff.close();
}

/**
* Write out an Eigen derived class to a file fn, with a boost::format string fmt.
* This assumes a file format with ii, jj, M(ii, jj)
* Assumes a symmetry condition applies, so the col loop only runs over
*   irow --> ncol
*/
template <typename DerivedBase>
void writeSymmMatrix_withIndices(std::string fn, std::string fmt, const Eigen::DenseBase<DerivedBase>& mat) {
    double nrow, ncol;
    nrow = mat.rows();
    ncol = mat.cols();

    // Open file
    std::ofstream ff(fn.c_str());

    if (ff.is_open()) {
        for (int ii=0; ii < nrow; ++ii)
            for (int jj = ii; jj < ncol; ++jj)
                ff << boost::format(fmt) % ii % jj % mat(ii, jj);
    } else {
        throw "Unable to open file\n";
    }

    // Close file
    ff.close();
}


/**
* Write out an Eigen derived class to a file fn, with a boost::format string fmt.
* This inserts an endline after each row.
*/
template <typename DerivedBase>
void writeMatrix(std::string fn, std::string fmt, const Eigen::DenseBase<DerivedBase>& mat) {
    double nrow, ncol;
    nrow = mat.rows();
    ncol = mat.cols();

    // Open file
    std::ofstream ff(fn.c_str());

    if (ff.is_open()) {
        for (int ii=0; ii < nrow; ++ii) {
            for (int jj = 0; jj < ncol; ++jj)
                ff << boost::format(fmt) % mat(ii, jj);
                ff << std::endl;
        }
    } else {
        throw "Unable to open file\n";
    }

    // Close file
    ff.close();
}


#endif // EIGEN_UTILS_H
