#ifndef GSLHIST_H
#define GSLHIST_H

#include <gsl/gsl_histogram.h>
#include "multi_array_wrap.h"

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

#endif // GSLHIST_H
