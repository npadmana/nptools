/* Wrappers around the Seo-Eisenstein code for computing distance errors
  We build this around the DETF code, not the full FoMSWG cosmological model,
  just as with the RSD code.
*/

#ifndef BAO_FORECAST_H
#define BAO_FORECAST_H

#include <cosmoutils.h>
#include <Eigen/Dense>
#include "wfirst_detf.h"

// Define the different Sigma's. Note that these are defined for the
// basic DETF cosmology, in particular, sigma8 = 0.8. So use accordingly, and
// be prepared for silly answers if not.
double Sigma_perp(double a);
double Sigma_par(double a);

// Basic BAO forecast code
Eigen::Matrix2d bao_forecast (
        const double number_density,   /* The number density in h^3 Mpc^-3 */
        const double sigma8,           /* The real-space, linear clustering amplitude */
        const double Sigma_perp,       /* The transverse rms Lagrangian displacement */
        const double Sigma_par,        /* The line of sight rms Lagrangian displacement */
        const double Sigma_z,          /* The line of sight rms comoving distance error due to redshift uncertainties */
                /* Note that Sigma_perp and Sigma_par are for pairwise differences,
                   while Sigma_z is for each individual object */
        const double beta,             /* The redshift distortion parameter */
        const double volume            /* The survey volume in h^-3 Gpc^3, set to 1 if input <=0 */
    );


// Overloaded for the more standard cases, where we specify a zmin and a zmax
Eigen::Matrix2d bao_forecast_shell (
        const double number_density,   /* The number density in h^3 Mpc^-3 */
        const double sigma8,           /* The real-space, linear clustering amplitude */
        const double Sigma_z,          /* The line of sight rms comoving distance error due to redshift uncertainties */
                /* Note that Sigma_perp and Sigma_par are for pairwise differences,
                   while Sigma_z is for each individual object */
        const double beta,             /* The redshift distortion parameter */
        const double zmin, const double zmax, const double area, /* zmin, zmax, area in deg^2 */
        const double recon=1.0,
        const bool isnumber = false
    );

#endif // BAO_FORECAST_H
