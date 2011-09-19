#include <cosmoutils.h>
#include <iostream>
#include <boost/format.hpp>
#include "wfirst_detf.h"

using namespace std;

int main(int argc, char** argv) {

    detf cosmo = detf_fiducial();

    // Print out useful bits of information
    double h = hubble(1.0, cosmo);
    cout << boost::format("h = %1$8.4f\n") % h;
    cout << boost::format("Omega_M = %1$7.4f\nOmega_DE = %2$7.4f\n") % (cosmo[4]/(h*h)) % cosmo[2];

    double aa;
    for (double zz = 0.1; zz < 3.0; zz += 0.1) {
        aa = 1./(1+zz);
        cout << boost::format("%1$4.2f %2$8.2f %3$8.2f %4$8.2f %5$8.4f\n")
                  % zz
                  % (comdis(aa, cosmo)*3000.0)
                  % (lumdis(aa, cosmo)*3000.0)
                  % (angdis(aa, cosmo)*3000.0)
                  % dm(aa, cosmo);
    }

}

