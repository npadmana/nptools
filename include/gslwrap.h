/** \file
 * Wrappers for GSL */

#ifndef GSLWRAP_H_
#define GSLWRAP_H_ 

#include <cmath>
// Boost functional includes
#include <boost/function.hpp>

// GSL includes
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_deriv.h>

#include "multi_array_wrap.h"

// Useful typedefs
typedef boost::function<double(double)> func_d_d;


// Adaptors for gsl_function
template<class F>
static double gslFunctionAdapter( double x, void* p)
{
    // Here I do recover the "right" pointer, safer to use static_cast
    // than reinterpret_cast.
    F* function = static_cast<F*>( p );
    return (*function)( x );
}

template<class F>
gsl_function convertToGslFunction( const F& f )
{
    gsl_function gslFunction;

    const void* p = &f;
    gslFunction.function = &gslFunctionAdapter<F>;
    // Just to eliminate the const.
    gslFunction.params = const_cast<void*>(p); 

    return gslFunction;
}


template<class F, int N=1000, int abs=-50, int rel=-20>
class Integrate {
    private :
        gsl_integration_workspace * wk;
        gsl_function ff;
        double epsabs, epsrel;
    public :
      Integrate(const F& f) {
          ff = convertToGslFunction(f);
          wk = gsl_integration_workspace_alloc(N);
          epsabs = exp(abs); epsrel = exp(rel);
      }
      ~Integrate() {
          gsl_integration_workspace_free(wk);
      }
      double operator()(double lo, double hi) {
          double result, err;
          gsl_integration_qags(&ff, lo, hi, epsabs, epsrel, N, wk, &result, &err);
          return result;
      }
};


/* Derivatives  - dir specifies the type, 0 = central, 1=forward, -1=reverse */
template<class F, int dir=0>
class Diff {
    private :
        gsl_function ff;
    public :
      Diff(const F& f) {
          ff = convertToGslFunction(f);
      }
      double operator()(double x, double h) {
          double result, err;
          switch (dir) {
            case 0 :
              gsl_deriv_central(&ff, x, h, &result, &err);
              break;
            case 1 :
              gsl_deriv_forward(&ff, x, h, &result, &err);
              break;
            case -1 :
              gsl_deriv_backward(&ff, x, h, &result, &err);
            }
          return result;
      }
};


class GSLHist
{
private :
       boost::shared_ptr<gsl_histogram> chist, whist;
       void alloc(int n);
public:
// Just use the default constructor
// Setting the range, will actually allocate the data
// This allows for easy subclassing of this object

       // Uniform bins
       void setbins(double xmin, double xmax, int n);
       // Arbitrary bins
       void setbins(const Eigen::VectorXd& bins);


       // Add to histogram
       void add(double x, double weight=1.0);

       // Get values
       Eigen::MatrixXd val();

       // reset
       void reset();
};


#endif /* GSLWRAP_H_ */
