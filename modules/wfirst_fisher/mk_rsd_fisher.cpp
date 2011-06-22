#include <iostream>
#include <string>
#include <cmath>
#include <tclap/CmdLine.h>
#include <eigen3/Eigen/Dense>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <cosmoutils.h>
#include <npio.h>
#include "wfirst_detf.h"
#include "zspace.h"
// We will need some more

using namespace std;
using namespace Eigen;

MatrixXd mk_lnfD_fisher(double aa, double ivar, detf cosmo0) {
    MatrixXd fish(ndetf, ndetf);
    vector<double> dtable = mk_deriv_table(lnfD, aa, ndetf, cosmo0);

    for (int ii=0; ii < ndetf; ++ii)
        for (int jj=ii; jj < ndetf; ++jj)
            fish(ii, jj) = fish(jj, ii) = ivar*dtable[ii]*dtable[jj];

     return fish;
}

int main(int argc, char** argv) {

   // Read in the command line
   TCLAP::CmdLine cmd("Generate an RSD fisher matrix",' ', "0.0");
   TCLAP::ValueArg<double> fsky("f", "fsky", "fraction of sky", true, 1.0, "fraction of sky [0.0, 1.0]");
   TCLAP::ValueArg<double> sigma8("s", "sigma8", "sigma8", false, 0.8, "sigma8");
   TCLAP::UnlabeledValueArg<string> infn("infn", "RSD input file (z, dz, nbar, bias)", true, "","RSD input filename (z, dz, nbar, bias)");
   TCLAP::UnlabeledValueArg<string> outfn("outfn", "RSD fisher file", true, "","RSD fisher filename");
   cmd.add(fsky);
   cmd.add(infn);
   cmd.add(outfn);
   cmd.parse(argc, argv);

   // Read in the file
   typedef boost::tuple<double, double, double, double> myrec;
   vector<myrec> ll = readAsciiFile(infn.getValue(), &TupleAdaptor<myrec, 4>);
   cout << ll.size() << " lines read in ....." << endl;
   double s8 = sigma8.getValue();
   cout << "Using sigma8 = " << s8 << endl;
   cout << "Using fsky = " << fsky.getValue() << endl;
   const double pi = 4.0*atan(1.0);
   double area = fsky.getValue() * (4.0*pi) * (180.0/pi) * (180.0/pi);
   cout << ".. corresponding to an area (in sq. deg.) = " << area << endl;

   // Define the fisher matrix
   MatrixXd dfish(ndetf, ndetf);
   dfish.setZero();
   detf fid = fiducial();


   // Now loop over all the entries, generating values
   double aa, amin, amax, z0, dz, nbar1, bias1, fval, vol, ivar;
   MatrixXd cov(2, 2); VectorXd nvec(1), bvec(1);
   BOOST_FOREACH( myrec l1, ll) {
        z0 = l1.get<0>(); dz = l1.get<1>(); nbar1 = l1.get<2>(); bias1 = l1.get<3>();
        aa = z2a(z0); amin = z2a(z0 - dz/2.0); amax = z2a(z0 + dz/2.0);
        vol = shellVol_Mpc_h(amin, amax, area, fid);
        fval = fgrowth(aa, fid);
        nvec[0] = nbar1; bvec[0] = bias1;
        zspace_mbias_pk(nvec, s8, bvec, fval, 0.0, vol, 0.1, cov);
        ivar = pow(fval*s8,2)/cov(1,1);
        cout << boost::format("%1$4.2f %2$4.2f %3$4.2f %4$6.3f %5$8.4e %6$8.4e\n") % z0 % dz % bias1 % fval % vol % ivar;
        dfish += mk_lnfD_fisher(aa, ivar, fid);
   }

   // Now marginalize the SN paramater
   dfish  = marginalizeSNparam(dfish);

   // Write this out
   writeDETFFisher(outfn.getValue(), dfish);

}
