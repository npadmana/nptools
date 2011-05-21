/* A simple example that reads in the 4 FOMSWG fisher matrices, and produces the DETF figure of merit.
 * Mostly as a check that I'm doing everything correctly.
 */


#include <iostream>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include "fisher_utils.h"
#include "fomswg.h"

using namespace std;
using namespace Eigen;

int main() {
    vector<string> fns;
    MatrixXd ff(nfswg, nfswg);
    ff.setZero();

    fns.push_back("PLANCK.dat");
    fns.push_back("preJDEM_SN.dat");
    fns.push_back("preJDEM_BAO.dat");
    fns.push_back("preJDEM_WL.dat");

    BOOST_FOREACH(string fn, fns) {
        ff += readFomSWG(fn);
    }

    // Remove the growth parameters
    ff = delete_parameter(ff, 5);
    ff = delete_parameter(ff, 6);

    // Set up marginalization
    vector<int> params;
    for (int ii = 0; ii<7; ++ii) params.push_back(ii);

    ff = marginalize(ff, params);
    Matrix2d detf;
    detf.setZero();
    double dwidwa, dwjdwa;
    for (int ii = 0; ii < 36; ++ii) {
        dwidwa = 0.025*ii +  0.0125;
        for (int jj = 0; jj < 36; ++jj) {
            dwjdwa = 0.025*jj + 0.0125;
            detf(0, 0) += ff(ii, jj);
            detf(1, 1) += ff(ii, jj) * dwjdwa * dwidwa;
            detf(0, 1) += ff(ii, jj) * dwjdwa;
            detf(1, 0) += ff(ii, jj) * dwidwa;
        }
    }


    cout << "The DETF Fisher matrix is :" << endl << detf << endl;
    cout << boost::format("The DETF FoM is %1%") % sqrt(detf.determinant()) << endl;



    cout << "Done...\n";
}
