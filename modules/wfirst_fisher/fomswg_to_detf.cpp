/* Convert the Stage III FoMSWG matrices to DETF */
#include <iostream>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include "fisher_utils.h"
#include "fomswg.h"
#include "wfirst_detf.h"

using namespace std;
using namespace Eigen;

int main() {
    vector<string> fns;
    MatrixXd ff1(nfswg, nfswg), ff(nfswg, nfswg);
    MatrixXd dfish1(ndetf, ndetf), dfish(ndetf-ndetf_nuis, ndetf-ndetf_nuis); // dfish accumulates only the non-SN parameters
    MatrixXd trans(nfswg, ndetf);
    ff.setZero();
    dfish.setZero();

    // Set up the list of files to read.
    fns.push_back("PLANCK");
    fns.push_back("preJDEM_SN");
    fns.push_back("preJDEM_BAO");
    fns.push_back("preJDEM_WL");

    // Transform matrix
    trans = mkTransformMatrix();


    BOOST_FOREACH(string fn, fns) {
        ff1 = readFomSWG(fn+".dat");
        ff += ff1;
        dfish1 = trans.transpose() * ff1 * trans;
        dfish1 = marginalizeSNparam(dfish1);
        writeDETFFisher(fn+"_detf.dat", dfish1);
        dfish += dfish1;
    }

    Matrix2d detf1, detf2;

    // First do this for the original FoMSWG matrices
    {
        // Remove the growth parameters
        ff = delete_parameter(ff, 5);
        ff = delete_parameter(ff, 6);
        // Set up marginalization
        vector<int> params;
        for (int ii = 0; ii<7; ++ii) params.push_back(ii);
        ff = marginalize(ff, params);
        detf1.setZero();
        double dwidwa, dwjdwa;
        for (int ii = 0; ii < 36; ++ii) {
            dwidwa = 0.025*ii +  0.0125;
            for (int jj = 0; jj < 36; ++jj) {
                dwjdwa = 0.025*jj + 0.0125;
                detf1(0, 0) += ff(ii, jj);
                detf1(1, 1) += ff(ii, jj) * dwjdwa * dwidwa;
                detf1(0, 1) += ff(ii, jj) * dwjdwa;
                detf1(1, 0) += ff(ii, jj) * dwidwa;
            }
        }
    }


    // Now do this for the new DETF Fisher matrices
    {
        // Remove the growth parameters
        dfish = delete_parameter(dfish, 8);
        // set up marginalization
        vector<int> params;
        for (int ii = 2; ii < 8; ++ii) params.push_back(ii);
        detf2 = marginalize(dfish, params);
    }


    cout << "The DETF Fisher matrix from the FOMSWG matrices is :" << endl << detf1 << endl;
    cout << "The DETF Fisher matrix from the DETF matrices is :" << endl << detf2 << endl;
    cout << "The DETF covariance matrix from the FOMSWG matrices is :" << endl << detf1.inverse() << endl;
    cout << "The DETF covariance matrix from the DETF matrices is :" << endl << detf2.inverse() << endl;
    cout << boost::format("The DETF FoM from FOMSWG is %1%") % sqrt(detf1.determinant()) << endl;
    cout << boost::format("The DETF FoM from DETF is %1%") % sqrt(detf2.determinant()) << endl;




    cout << "Done...\n";
}
