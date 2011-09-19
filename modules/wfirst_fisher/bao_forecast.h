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




#endif // BAO_FORECAST_H
