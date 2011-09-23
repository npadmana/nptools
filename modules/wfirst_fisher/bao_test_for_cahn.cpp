// Some code for tests that Bob Cahn would like for the BAO data
// Nikhil Padmanabhan, Yale
// Sitting at LBL, Sept 2011

#include <iostream>
#include <string>
#include <cmath>
#include <gslwrap.h>
#include <tclap/CmdLine.h>
#include <eigen3/Eigen/Dense>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <cosmoutils.h>
#include <npio.h>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include "wfirst_detf.h"
#include "bao_forecast.h"
// We will need some more

using namespace std;
using namespace Eigen;
using namespace boost::lambda;

double lnHs(double a, detf c) {
    return log(hubble(a, c) * sound_horizon_eh98_fit(c));
}

double lnDa_s(double a, detf c) {
    return log(propmotdis(a, c) / sound_horizon_eh98_fit(c));
}


double _mkderiv_impl(double x, int idetf, int iobs, double amid, detf c) {
   c[idetf] = x;
   if (iobs == 0) {
       return lnDa_s(amid, c);
   } else {
       return lnHs(amid, c);
   }
}

MatrixXd MkDerivs(double amid, detf c0) {
    MatrixXd deriv(ndetf, 2);

    for (int idetf=0; idetf < ndetf; ++idetf)
        for (int iobs=0; iobs < 2; ++iobs) {
            deriv(idetf, iobs) = Diff<func_d_d>( bind(_mkderiv_impl, _1, idetf, iobs, amid, c0))(c0[idetf], 1.e-5);
        }

    return deriv;
}



int main(int argc, char** argv) {

   // Read in the command line
   TCLAP::CmdLine cmd("Generate an RSD fisher matrix",' ', "0.0");
   TCLAP::ValueArg<double> area("a", "area", "area of sky [deg^2]", true, 1.e4, "area of sky");
   TCLAP::ValueArg<double> sigma8("s", "sigma8", "sigma8 ", false, 0.8, "sigma8 [0.8]");
   TCLAP::ValueArg<double> bias0("b", "bias", "bias at z=0", false, 1.0, "bias [1.0]");
   TCLAP::ValueArg<double> dz("z", "dz", "delta z [0.1]", false, 0.1, "delta z");
   TCLAP::UnlabeledValueArg<string> infn("infn", "BAO input file (zmin, zmax, dn/dzdA)", true, "","BAO input file (zmin, zmax, dn/dzdA)");
   TCLAP::UnlabeledValueArg<string> outfn("outfn", "BAO fisher file", true, "","BAO fisher filename");
   cmd.add(area);
   cmd.add(sigma8);
   cmd.add(bias0);
   cmd.add(dz);
   cmd.add(infn);
   cmd.add(outfn);
   cmd.parse(argc, argv);

   // Define the basic cosmology
   detf fidcosmo = detf_fiducial();
   double D0 = growth(1.0, fidcosmo);


   // Read in the file
   typedef boost::tuple<double, double> myrec;
   vector<myrec> ll = readAsciiFile(infn.getValue(), &TupleAdaptor<myrec, 2>);
   cout << ll.size() << " lines read in ....." << endl;
   double s8 = sigma8.getValue();
   cout << "Using sigma8 = " << s8 << endl;
   cout << "Using area = " << area.getValue() << endl;
   cout << "Using dz = " << dz.getValue() << endl;
   cout << "Using bias(z=0) = " << bias0.getValue() << endl;

   double zmin, zmax, zmid, num, errD, errH, bb, beta, amid, s8z, Dz;
   Matrix2d cov;
   MatrixXd dmat(ndetf, 2);
   MatrixXd dfish(ndetf, ndetf);
   dfish.setZero();
   cout << "#zmid num bias beta errD(%) errH(%) \n";
   BOOST_FOREACH( myrec l1, ll) {
        zmid = l1.get<0>(); num = l1.get<1>(); amid = z2a(zmid);
        Dz = growth(amid, fidcosmo)/D0;
        s8z = s8 * Dz; // Update sigma8
        bb = bias0.getValue()/Dz; // Update bias
        beta = fgrowth(amid, fidcosmo)/bb;
        zmin = zmid - dz.getValue()/2.0; zmax = zmid + dz.getValue()/2.0;
        cov = bao_forecast_shell(num, bb, s8z, 0.0, beta, zmin, zmax, area.getValue(), 0.5, true);
        errD = sqrt(cov(0,0)); errH = sqrt(cov(1,1));
        cov = cov/1.e4;
        cout << boost::format("%1$4.2f %2$7.1f %3$4.2f %4$5.3f %5$7.3f %6$7.3f\n") % zmid % num % bb % beta % errD % errH;
        dmat =  MkDerivs(amid, fidcosmo);
        dfish += dmat * (cov.inverse() * dmat.transpose());
   }

   // Now marginalize the SN paramater
   dfish  = marginalizeSNparam(dfish);

   // Write this out
   writeDETFFisher(outfn.getValue(), dfish);

}
