#include <iostream>
#include <cmath>
#include <tclap/CmdLine.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <eigen3/Eigen/Dense>
#include "fisher_utils.h"
#include "wfirst_detf.h"


using namespace std;
using namespace boost;
using namespace Eigen;

int main(int argc, char** argv) {
   TCLAP::CmdLine cmd("Combine DETF fisher matrices",' ', "0.0");
   TCLAP::UnlabeledMultiArg<string> fns("fns", "fisher matrix files", true, "filenames");
   cmd.add(fns);
   cmd.parse(argc, argv);

   MatrixXd dfish(ndetf-ndetf_nuis, ndetf-ndetf_nuis), dtemp;
   dfish.setZero();

   BOOST_FOREACH (string s, fns.getValue()) {
       cout << boost::format("Processing file.... %1%\n") % s;
       dfish += readDETF(s);
   }

   double FoM, FoMwg, FoMg;
   // Compute DETF FoM
   {
    vector<int> params;
    dtemp = delete_parameter(dfish, 8);
    for (int ii=2; ii<8; ++ii) params.push_back(ii);
    dtemp = marginalize(dtemp, params);
    FoM = sqrt(dtemp.determinant());
   }

   // Allow gamma to vary
   {
    vector<int> params;
    for (int ii=2; ii<9; ++ii) params.push_back(ii);
    dtemp = marginalize(dfish, params);
    FoMwg = sqrt(dtemp.determinant());
   }

   // Allow gamma FoM
   {
    vector<int> params;
    for (int ii=0; ii<8; ++ii) params.push_back(ii);
    dtemp = marginalize(dfish, params);
    FoMg = dtemp(0,0);
   }

   cout << boost::format("The DETF FoM (gamma fixed) = %1$6.2f \n") % FoM;
   cout << boost::format("The DETF FoM (gamma varying) = %1$6.2f \n") % FoMwg;
   cout << boost::format("The gamma FoM  = %1$6.2f \n") % FoMg;

}
