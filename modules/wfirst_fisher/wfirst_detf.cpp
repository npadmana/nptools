#include "wfirst_detf.h"
#include <cmath>
#include <gslwrap.h>
#include <npio.h>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include "fisher_utils.h"

using namespace std;
using namespace Eigen;
using namespace boost::lambda;

detf fiducial() {
    detf x;
    x[0] = -1.0; //w0
    x[1] = 0.0; //wa
    x[2] = 0.7436; // Omega_DE
    x[3] = 0.0; // Omega_K
    x[4] = 0.1326; // omega_mh^2
    x[5] = 0.0227; // omega_bh^2
    x[6] = 0.963; // ns
    x[7] = -19.9628; // ln Delta
    x[8] = 0.0; // Delta gamma
    x[9] = 0.0; // Script M

    return x;
}

double Omega_K(detf x) {return x[3];}

double hubble(double a, detf x) {
    double OmegaM0 = 1.0 - x[2] - x[3];
    double h20 = x[4]/OmegaM0;

    double hh;

    // Omega M, and Omega K
    hh = x[4]/pow(a, 3) + x[3]*h20/pow(a, 2);

    // DE
    double val = -3.0 * ((1.0+x[0] + x[1]) * log(a) - (1-a)*x[1]);
    hh += x[2]*h20*exp(val);

    return hh;
}

Fswg detf2Fswg (detf c) {
    Fswg x;
    double OmegaM0 = 1.0 - c[2] - c[3];
    double h20 = c[4]/OmegaM0;

    x[0] = c[6]; //ns
    x[1] = c[4]; //omh2
    x[2] = c[5]; //ombh2
    x[3] = c[3]*h20; // K
    x[4] = c[2]*h20; // DE
    x[5] = c[8]; // gamma
    x[6] = c[9]; // M
    x[7] = 0.0; // G is fixed
    x[8] = c[7]; // ln D

    // Now we do the DE EoS parameters
    double aa;
    for (int ii=9; ii < nfswg; ++ii) {
        aa = 1.0 - 0.025*ii - 0.0125;
        x[ii] = c[0] + (1.0-aa)*c[1];
    }

    return x;
}


double _detf2Fswg(double x, detf c, int ii, int jj) {
    c[ii] = x;
    return detf2Fswg(c)[jj];
}

MatrixXd mkTransformMatrix() {
    detf cc = fiducial();

    MatrixXd transform(nfswg, ndetf);

    for (int ii = 0; ii < nfswg; ++ii)
        for (int jj = 0; jj < ndetf; ++jj)
            transform(ii, jj) = Diff<func_d_d>( bind(_detf2Fswg, _1, cc, jj, ii))(cc[jj], 1.e-5);

    return transform;
}


MatrixXd readDETF(string fn) {
    typedef boost::array<double, 9> rec;
    vector<rec> ll;

    ll = readAsciiFile(fn, &ArrayAdaptor<rec, 9>);

    MatrixXd fish(ndetf-ndetf_nuis, ndetf-ndetf_nuis);

    for (int ii=0; ii < ndetf-ndetf_nuis; ++ii)
        for (int jj=0; jj< ndetf-ndetf_nuis; ++jj)
            fish(ii, jj) = ll[ii][jj];

    return fish;
}


void writeDETFFisher(std::string fn, const MatrixXd& mat) {
    writeMatrix(fn, " %1$25.15e", mat);
}



Eigen::MatrixXd marginalizeSNparam(const Eigen::MatrixXd& mat) {
    vector<int> params;
    params.push_back(9);
    return marginalize(mat, params);
}


double snmag(double a, detf cosmo) {return dm(a, cosmo) + cosmo[9];}

