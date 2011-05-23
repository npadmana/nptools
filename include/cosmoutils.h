/** Templated cosmological functions. These can use any baseline cosmology, as long as they
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
  double OmK = Omega_K(cosmo);
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
template <class C> double dm(double a, C cosmo) {return 5.0 * log10(lumdis(a, cosmo) + 25.0);}


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
