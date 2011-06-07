#ifndef FORECAST_H
#define FORECAST_H

#include <cosmoutils.h>
#include <Eigen/Dense>
#include "fomswg.h" // We will necessarily need to modify fomswg
#include "bao_forecast.h"

// These are defined for \sigma_8 = 0.8! If you use a different value, be prepared for issues!
double SigmaPerp(double a);
double SigmaPar(double a);


class BAO_Forecast_ZSlice {
    public :
        // Store everything here
        double zmin, zmax, area, vol;
        double nbar, bias, beta;
        double sperp, sll, sigmaz;
        double recon;
        Matrix2d cov;

        BAO_Forecast_ZSlice(double _zmin, double _zmax, double _area, double _nbar, double _bias, double _recon=1.0, double _sigmaz=0.0);

        // Helper function that calculates all the derived quantities
        void derived();

};



#endif // FORECAST_H
