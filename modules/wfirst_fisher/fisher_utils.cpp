#include "fisher_utils.h"
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <iostream>

using namespace std;
using namespace Eigen;

MatrixXd delete_parameter(const MatrixXd& in, int par) {
    return delete_row(delete_col(in, par), par);
}


MatrixXd delete_row(const MatrixXd& in, int irow) {
    int rows = in.rows();
    int cols = in.cols();

    MatrixXd out(rows-1, cols);

    // Extract
    out.block(0, 0, irow, cols) = in.block(0, 0, irow, cols);
    out.block(irow, 0, rows-irow-1, cols) = in.block(irow+1, 0, rows-irow-1, cols);

    return out;
}

MatrixXd delete_col(const MatrixXd& in, int icol) {
    int rows = in.rows();
    int cols = in.cols();

    MatrixXd out(rows,  cols-1);

    // Extract
    out.block(0, 0, rows, icol) = in.block(0, 0, rows, icol);
    out.block(0, icol, rows, cols-icol-1) = in.block(0, icol+1, rows, cols-icol-1);

    return out;
}

// params is passed by value, since we will modify it.
MatrixXd marginalize(const MatrixXd& in, vector<int> params) {
    // Set up the sets
    //    full -- everything
    //    params -- the nuisance parameters
    //    rem  -- what remains
    int npar= in.rows();
    int nnuis = params.size();
    int nrem = npar - nnuis;
    vector<int> full(npar), rem(nrem);
    vector<int>::iterator it;

    for (int ii=0; ii< npar; ++ii) full[ii] = ii;
    sort(params.begin(), params.end());
    it = set_difference(full.begin(), full.end(), params.begin(), params.end(), rem.begin());

    if (nrem != (it - rem.begin()))cout << "Incorrect number of parameters obtained!\n";
    cout << boost::format("Initial = %1%, Marginalized = %2%, Remaining %3% \n") % npar % nnuis %  nrem;

    // Now set up the matrices that we need
    MatrixXd Frr, Fqq, Fqr; // This is the FOMSWG notation
    Frr = Fqq = Fqr = in;

    // Build
    for (int ii = 0; ii< nnuis; ++ii) {
        Fqq = delete_parameter(Fqq, params[ii] - ii); // the -ii is to correctly handle the fact that all the columns and rows are decreasing
        Fqr = delete_row(Fqr, params[ii] -ii);
    }

    for (int ii = 0; ii < nrem; ++ii) {
        Frr = delete_parameter(Frr, rem[ii] - ii);
        Fqr = delete_col(Fqr, rem[ii] - ii);
    }

    // Diagonalize Frr
    SelfAdjointEigenSolver<MatrixXd> eigen(Frr);
    VectorXd dd = eigen.eigenvalues();
    Fqr = Fqr*eigen.eigenvectors();
    cout << Fqr.rows() << " " << Fqr.cols() << endl;
    for (int ii=0; ii < nnuis; ++ii) {
        if (dd[ii] < 1.e-5) {
            cout << boost::format("WARNING! %1% parameter has a very small eigenvalue....\n")% ii;
            cout << "Printing out the corresponding terms for you to check.....\n";
            cout << Fqr.col(ii) << endl;
            dd[ii] = 0.0;
        } else {
            dd[ii] = 1./dd[ii];
        }
    }

    Fqq = Fqq - Fqr * dd.asDiagonal() * Fqr.transpose();

    return Fqq;
}


