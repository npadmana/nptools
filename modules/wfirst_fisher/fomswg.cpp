#include "fomswg.h"
#include "npio.h"

using namespace std;
using namespace Eigen;


// Generate the fiducial cosmology
Fswg Fiducial() {
  Fswg arr;

  // Do this by hand
  arr[0] = 0.963; //ns
  arr[1] = 0.1326; // omega_m
  arr[2] = 0.0227; // omega_b
  arr[3] = 0; // omega_k
  arr[4] = 0.3844; // omega_de
  arr[5] = 0.0; // Delta gamma
  arr[6] = 0.0; // Delta M
  arr[7] = 0.0; // ln G0
  arr[8] = -19.9628; // ln \Delta_{\zeta}^2

  // Now set the DE EoS parameters
  for (int ii=9; ii < nfswg; ++ii) arr[ii] = -1.0;

  return arr;
}

double Omega_K(Fswg cosmo) {return cosmo[3]/(cosmo[1] + cosmo[3] + cosmo[4]);}

// Return the hubble parameter at scale factor a
// in units of 100 km/s/Mpc
double hubble(double a, Fswg cosmo) {
  double hh;

  hh = cosmo[1]/pow(a, 3) + cosmo[3]/(a*a); // No radiation

  // Now handle the DE term.
  double astart = 1.0, aend = 1.0, eosexp = 0.0, alo;
  int ii = startw;
  while ((a < aend) && (ii < nfswg)) {
    astart = 1.0 - 0.025*(ii-startw);
    aend = astart - 0.025;
    alo = (aend < a) ? a : aend; // Set up the lower integral bound
    eosexp += (log(astart) - log(alo))*3.0*(1.0 + cosmo[ii]);
    ii++;
  }
  hh += cosmo[4] * exp(eosexp);

  return sqrt(hh);
}

double snmag(double a, Fswg cosmo) {return dm(a, cosmo) + cosmo[6];}


MatrixXd readFomSWG(string fn) {
    typedef boost::tuple<int, int, double> rec;
    vector<rec> ll;

    ll = readAsciiFile(fn, &TupleAdaptor<rec, 3>);

    MatrixXd fish(nfswg, nfswg);
    fish.setZero();

    BOOST_FOREACH (rec ii, ll) {
        fish(ii.get<0>(), ii.get<1>()) = ii.get<2>();
        fish(ii.get<1>(), ii.get<0>()) = ii.get<2>();
    }

    return fish;
}
