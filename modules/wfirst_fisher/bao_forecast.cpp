#include "bao_forecast.h"

// Define the fiducial cosmology
detf _fidcosmo = detf_fiducial();

double Sigma_perp(double a) {
    return 1.0;
}

double Sigma_par(double a) {
    return 1.0;
}
