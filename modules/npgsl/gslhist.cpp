#include "gslwrap.h"

using namespace std;
using namespace Eigen;

void GSLHist::alloc(int n) {
    chist = boost::shared_ptr<gsl_histogram>( gsl_histogram_alloc(n), std::ptr_fun(gsl_histogram_free));
    whist = boost::shared_ptr<gsl_histogram>( gsl_histogram_alloc(n), std::ptr_fun(gsl_histogram_free));
}


void GSLHist::setbins(double xmin, double xmax, int n) {
    alloc(n);
    gsl_histogram_set_ranges_uniform(chist.get(), xmin, xmax);
}

void GSLHist::setbins(const VectorXd& bins) {
    alloc(bins.size());
    gsl_histogram_set_ranges(chist.get(), bins.data(),  bins.size());
}


void GSLHist::add(double x, double weight) {
    gsl_histogram_accumulate(whist.get(), x, weight);
    gsl_histogram_increment(chist.get(), x);
}


MatrixXd GSLHist::val() {
   int nbins = gsl_histogram_bins(chist.get());
   MatrixXd retval(nbins, 4);

   // Fill 'er in
   for (int ii=0; ii< nbins; ++ii) {
       retval(ii, 2) = gsl_histogram_get(whist.get(), ii);
       retval(ii, 3) = gsl_histogram_get(chist.get(), ii);
       gsl_histogram_get_range(whist.get(), ii, &retval(ii,0), &retval(ii, 1));
   }
   return retval;
}

void GSLHist::reset() {
    gsl_histogram_reset(chist.get());
    gsl_histogram_reset(whist.get());
}
