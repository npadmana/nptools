#include <iostream>
#include <vector>
#include <string>
#include <tclap/CmdLine.h>
#include <eigen3/Eigen/Dense>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <cosmoutils.h>
#include <npio.h>
#include "wfirst_detf.h"
// We will need some more

using namespace std;
using namespace Eigen;

// Simple fisher matrix generator for SN
MatrixXd mk_sn_fisher(double aa, double ivar, detf cosmo0) {
    MatrixXd fish(ndetf, ndetf);
    vector<double> dtable = mk_deriv_table(snmag, aa, ndetf, cosmo0);

    for (int ii=0; ii < ndetf; ++ii)
        for (int jj=ii; jj < ndetf; ++jj)
            fish(ii, jj) = fish(jj, ii) = ivar*dtable[ii]*dtable[jj];

     return fish;
}

int main(int argc, char** argv) {

   // Read in the command line
   TCLAP::CmdLine cmd("Generate a SN fisher matrix",' ', "0.0");
   TCLAP::UnlabeledValueArg<string> infn("infn", "SN error file", true, "","SN error filename");
   TCLAP::UnlabeledValueArg<string> outfn("outfn", "SN fisher file", true, "","SN fisher filename");
   cmd.add(infn);
   cmd.add(outfn);
   cmd.parse(argc, argv);

   // Read in the file
   typedef boost::tuple<double, double> myrec;
   vector<myrec> ll = readAsciiFile(infn.getValue(), &TupleAdaptor<myrec, 2>);

   // Define the fisher matrix
   MatrixXd dfish(ndetf, ndetf), dfish1(ndetf, ndetf);
   detf fid = fiducial();
   BOOST_FOREACH( myrec l1, ll) {
        double aa = 1./(1.+l1.get<0>());
        double ivar = 1./pow(l1.get<1>(),2);
        dfish1 = mk_sn_fisher(aa, ivar, fid);
        dfish += dfish1;
   }

   // Now marginalize the SN paramater
   dfish  = marginalizeSNparam(dfish);

   // Write this out
   writeDETFFisher(outfn.getValue(), dfish);

}
