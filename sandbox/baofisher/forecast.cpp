#include "forecast.h"

Fswg _fiducial = Fiducial(); // Set up a fiducial cosmology
static const double sigma8 = 0.8;
static const double _sigma0 = 12.4 * 0.8/0.9; // Define this

double SigmaPerp(double a) {
    return _sigma0 * growth(a, _fiducial);
}

double SigmaPar(double a) {
    return SigmaPerp(a)*(1+ fgrowth(a, _fiducial));
}

BAO_Forecast_ZSlice::BAO_Forecast_ZSlice(double _zmin, double _zmax, double _area, double _nbar, double _bias, double _recon, double _sigmaz) {
    // Basic params
    zmin = _zmin;
    zmax = _zmax;
    area = _area;

    // Galaxy params
    nbar = _nbar;
    bias = _bias;
    recon = _recon;
    sigmaz = _sigmaz;
    derived();
}

void BAO_Forecast_ZSlice::derived() {
    double zmid = (zmin + zmax)/2.0;
    double amid = z2a(zmid);
    vol = shellVol_Gpc_h(z2a(zmin), z2a(zmax), area, _fiducial);
    beta = fgrowth(amid, _fiducial)/bias;
    sperp = SigmaPerp(amid)*recon;
    sll = SigmaPar(amid)*recon;

    // Now do the actual call
    double Drms, Hrms, r, Rrms;

    bao_forecast(nbar, sigma8*bias, sperp, sll, sigmaz, beta, vol, &Drms, &Hrms, &r, &Rrms);
    cov << Drms*Drms << Drms * Hrms *r << Drms * Hrms *r  << Hrms * Hrms;
}


