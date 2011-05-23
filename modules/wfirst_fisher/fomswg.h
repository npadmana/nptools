/* This is an example of how you might define a cosmology.
 *
 * For the Figure of Merit Science Working Group cosmological definition.
 * Note that isn't abstract, but can be easily fed into any routine.
 */

#ifndef FOMSWG_H
#define FOMSWG_H

#include <string>
#include <boost/array.hpp>
#include <eigen3/Eigen/Dense>
#include <cosmoutils.h>
#include <npio.h>

const int nfswg = 45; // Number of FOMSWG parameters
const int startw = 9; // The starting index for the DE EoS parameters

// Note that this automagically means that Fswg has the [ ] operator
typedef boost::array<double, nfswg> Fswg;

// Generate the fiducial cosmology
Fswg Fiducial();

// Return Omega_K at z=0
double Omega_K(Fswg cosmo);

// Return the hubble parameter at scale factor a
// in units of 100 km/s/Mpc
// Note that this is an actual hubble constant, not in terms of h
double hubble(double a, Fswg cosmo);

/////////////////////////////////////////////////////////
// EVERYTHING ABOVE DEFINES A MINIMAL INTERFACE
// ANYTHING ELSE GOES BELOW HERE
/////////////////////////////////////////////////////////

// SN magnitudes.
double snmag(double a, Fswg cosmo);

Eigen::MatrixXd readFomSWG(std::string fn);

#endif // FOMSWG_H
