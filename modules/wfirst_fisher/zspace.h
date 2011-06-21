#ifndef ZSPACE_H
#define ZSPACE_H
#include <gsl/gsl_linalg.h>
#include <Eigen/Dense>

void zspace_mbias_pk  (
                      const int NSAMP, // number of samples
                      double *nbar,    // The number density in h^3 Mpc^-3
                      double sigma8,   // The real-space, linear clustering amplitude
                      double *bias,    // The real-space, linear bias
                      double f,	       // f ~ Omega_m^(0.6)
                      double Sigma_z,  // z error translated into comoving distance
                      double vol_mpc,  // The survey volume in h^-3 Mpc^3
                      double kmax,     // Maximum k-value to integrate to
                      Eigen::MatrixXd& invfish  // inverse covariance matrix for bs and fs
    );


#endif // ZSPACE_H
