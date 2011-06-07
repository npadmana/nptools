/** \file 
 * Templated cosmological functions. These can use any baseline cosmology, as long as they
 * follow a basic interface design.
 * See the example implementations for examples of how to do this.
 * This structure avoids the clumsiness of an fully object oriented approach, but keeps the
 * extensibility available.
 *
 * Nikhil Padmanabhan, Yale, May, 2011
 */

#ifndef COSMOUTILS_H
#define COSMOUTILS_H

#include <cmath>
#include <vector>
#include "gslwrap.h"

// Useful cosmological functions
inline double a2z(double a) {return 1./a - 1.0;}
inline double z2a(double z) {return 1./(1.0 +z);}


// Return Omega_M as a function of redshift
template <class C>
double Omega_M(double a, C cosmo) {
    double om0 = Omega_M0(cosmo); // Get the redshift zero measurement

    // scale by a**3
    om0 /= pow(a, 3);

    // The critical density is propto h^2
    double crit = hubble(a, cosmo)/hubble(1.0, cosmo);
    double crit2 = crit*crit;

    return om0/crit2;
}

// Return the logarithmic growth rate as a function of redshift.
// This assumes that a growth index gamma function exists for this cosmology
// This simply assume fgrowth = \Omega_M^gamma
template <class C>
double fgrowth(double a, C cosmo) {return pow(Omega_M(a, cosmo), grgamma(cosmo));}

// Helper class for growth function
template <class C>
class _growth {
    private :
        C cosmo;
    public :
        _growth(C _cosmo) : cosmo(_cosmo) {}
        double operator() (double a) {return (fgrowth(a, cosmo)-1)/a;}
};

// Compute the growth function
template <class C>
double growth(double a, C cosmo) {
  double tmp = Integrate<func_d_d> (_growth<C>(cosmo))(1.e-3, a);
  return exp(tmp)*a;
}



// Helper class for comdis
// This is necessary, because of the templating inside comdis, which prevents us from using lambdas
template <class C>
class _comdis {
    private :
        C cosmo;
    public :
        _comdis(C _cosmo) {cosmo = _cosmo;}
        double operator() (double a) {return 1./(a*a*hubble(a, cosmo));}
};

// Compute the comoving distance
template <class C>
double comdis(double a, C cosmo) {
  return Integrate<func_d_d> (_comdis<C>(cosmo))(a, 1.0); //Integrate da/a^2 H
}

// Proper motion distance
template <class C>
double propmotdis(double a, C cosmo) {
  double OmK = Omega_K0(cosmo);
  double hval = hubble(1.0, cosmo);
  double omK = OmK * hval * hval;
  double somK = sqrt(abs(omK));

  double r = comdis(a, cosmo); // Comoving distance.... in units of 3000 Mpc
                                  // Note, there is no h^-1 in this

  if ((somK * r) < 1.e-2)
       {return r*(1 + pow(somK*r, 2)/6.0 + pow(somK*r,4)/120.0);}
  else if (OmK < 0)
       {return sin(somK * r)/somK;}
  else {return sinh(somK * r)/somK;}
}

template <class C> double lumdis(double a, C cosmo) {return propmotdis(a, cosmo)/a;}
template <class C> double angdis(double a, C cosmo) {return propmotdis(a, cosmo)*a;}
template <class C> double dm(double a, C cosmo) {return 5.0 * log10(lumdis(a, cosmo)) + 25.0 + 5.0*log10(3000.0);}

// Return the spherical volume of a shell to amax
// See Hogg, Distance Measures in Cosmology
template <class C>
double vol(double a, C cosmo) {
    double OmK = Omega_K0(cosmo);
    double DM = comdis(a, cosmo)/hubble(1.0, cosmo);

    // Define some temporaries and compute 'em
    double f1, f2, sOmK, vol;
    sOmK = sqrt(abs(OmK));
    f1 = DM * sqrt(1.0 + OmK * DM * DM);
    f2 = sOmK * DM;

    // Now handle the various cases
    if (OmK < 1.e-5) {
        // Flat
        vol = 1./3.;
    } else if (OmK < 0) {
        vol = f1 - asin(f2)/sOmK;
    } else {
        vol = f1 - asinh(f2)/sOmK;
    }

    return vol;
}

// We define a whole suite of overloaded volume calculations

// This computes the volume of a shell with unit area in units of
// 3000^3 Gpc^3
template <class C>
double shellVol(double a0, double a1, C cosmo) {
    return vol(a1, cosmo) - vol(a0, cosmo);
}


// This computes the volume of a shell of area A in degrees^2
// in (Gpc/h)^3
template <class C>
double shellVol_Gpc_h(double a0, double a1, double area, C cosmo) {
    double str2deg2 = pow(45/atan(1.0), 2);
    area /= str2deg2;
    double h3 = pow(hubble(1.0, cosmo), 3); // h^3
    return area * 27.0* shellVol(a0, a1, cosmo)/h3;
}

// Same as previous, but in (Mpc/h)^3
template <class C>
double shellVol_Mpc_h(double a0, double a1, double area, C cosmo) {
    return shellVol_Gpc_h(a0, a1, area, cosmo) * 1.e9;
}














// Simple fisher matrix routines below.
// Using these requires a [ ] interface to the cosmological parameters

// Compute derivatives for the Fisher matrix
template <class C>
class _partialD {
    private :
       C cosmo;
       int i;
       double a;
       double (*observable) (double, C);
    public :
       _partialD(double (*_obs) (double, C), double _a, int _i, C _cosmo) { observable = _obs; cosmo = _cosmo; a = _a; i = _i;}
       double operator() (double x) {cosmo[i] = x; return observable(a, cosmo);}
};

template <class C>
std::vector<double> mk_deriv_table(double (*observable)(double, C), double aa, int Npar, C cosmo0) {
    std::vector<double> dtable(Npar);
    for (int ii=0; ii < Npar; ++ii) {
        _partialD<C> pd(observable, aa, ii, cosmo0);
        dtable[ii] = Diff < func_d_d > (pd) (cosmo0[ii], 1.e-5);
    }
    return dtable;
}



#endif // COSMOUTILS_H
