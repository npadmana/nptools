/* Convert the Stage III FoMSWG matrices to DETF */
#include <iostream>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <tclap/CmdLine.h>
#include "fisher_utils.h"
#include "fomswg.h"
#include "wfirst_detf.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
    TCLAP::CmdLine cmd("Convert FoMSWG fisher matrice to DETF format",' ', "0.0");
    TCLAP::UnlabeledMultiArg<string> fns("fns", "fisher matrix input and output files", true, "filenames");
    cmd.add(fns);
    cmd.parse(argc, argv);

    MatrixXd ff1(nfswg, nfswg), ff(nfswg, nfswg);
    MatrixXd dfish1(ndetf, ndetf), dfish(ndetf-ndetf_nuis, ndetf-ndetf_nuis); // dfish accumulates only the non-SN parameters
    MatrixXd trans(nfswg, ndetf);
    dfish.setZero();

    // Transform matrix
    trans = mkTransformMatrix();


    ff1 = readFomSWG(fns.getValue()[0]);
    dfish1 = trans.transpose() * ff1 * trans;
    dfish1 = marginalizeSNparam(dfish1);
    writeDETFFisher(fns.getValue()[1], dfish1);
    dfish += dfish1;


}
