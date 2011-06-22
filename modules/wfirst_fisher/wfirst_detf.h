/* The cosmology that the WFIRST SDT team has adopted.
 */
#ifndef WFIRST_DETF_H
#define WFIRST_DETF_H

#include <boost/array.hpp>
#include <eigen3/Eigen/Dense>
#include <cosmoutils.h>
#include <eigen_utils.h>
#include "fomswg.h"

const int ndetf = 10; // Include a script M for SN
const int ndetf_nuis = 1; // the SN param
typedef boost::array<double, ndetf> detf;

detf fiducial();
double Omega_K0(detf c);
double Omega_M0(detf c);
double grgamma(detf c);
double hubble(double a, detf c);

/*** The basic interface ends here ***/

double snmag(double a, detf cosmo);
double lnfD(double a, detf cosmo);
// Convert DETF to Fswg
Fswg detf2Fswg(detf c);
Eigen::MatrixXd mkTransformMatrix();

Eigen::MatrixXd marginalizeSNparam(const Eigen::MatrixXd& mat);

Eigen::MatrixXd readDETF(std::string fn);
void writeDETFFisher(std::string fn, const Eigen::MatrixXd& mat);

#endif // WFIRST_DETF_H
