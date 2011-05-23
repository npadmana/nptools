/* The cosmology that the WFIRST SDT team has adopted.
 */
#ifndef WFIRST_DETF_H
#define WFIRST_DETF_H

#include <boost/array.hpp>
#include <eigen3/Eigen/Dense>
#include <cosmoutils.h>
#include <eigen_utils.h>
#include <npio.h>
#include "fomswg.h"

const int ndetf = 10; // Include a script M for SN
typedef boost::array<double, ndetf> detf;

detf fiducial();
double Omega_K(detf c);
double hubble(double a, detf c);

/*** The basic interface ends here ***/

// Convert DETF to Fswg
Fswg detf2Fswg(detf c);
Eigen::MatrixXd mkTransformMatrix();

void writeDETFFisher(std::string fn, const Eigen::MatrixXd& mat);

#endif // WFIRST_DETF_H
