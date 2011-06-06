#include "linelist.h"
#include <stdexcept>
#include <boost/lambda/lambda.hpp>
#include <functional>
#include <algorithm>

/*
List of lines from the SDSS website -- not all may be implemented in the code
# wave   label
1033.30  OVI
1215.67  Ly_alpha
1239.42  NV
1305.53  OI
1335.52  CII
1399.8   SiIV+OIV
1545.86  CIV
1640.4   HeII
1665.85  OIII
1857.4  AlIII
1908.27  CIII
2326.0  CII
2439.5   NeIV
2800.32  MgII
3346.79  NeV
3426.85  NeV
3728.30  OII
3798.976 H_theta
3836.47  H_eta
3889.0   HeI
3934.777 K
3969.588 H
4072.3   SII
4102.89  H_delta
4305.61  G
4341.68  H_gamma
4364.436 OIII
4862.68  H_beta
4960.295 OIII
5008.240 OIII
5176.7   Mg
5895.6   Na
6302.046 OI
6365.536 OI
6549.86  NII
6564.61  H_alpha
6585.27  NII
6707.89  Li
6718.29  SII
6732.67  SII
*/


using namespace std;
using namespace Eigen;
using namespace boost::lambda;

LineList::LineList()
{
    // Insert lines here
    llist.insert(linemap::value_type("CIV", 1545.86));
    llist.insert(linemap::value_type("CIII", 1908.27));
    llist.insert(linemap::value_type("MgII", 2800.32));
}


int LineList::nlines(std::string linestr) {
    return llist.count(linestr);
}

double _get_second(pair<string, double> it) {
    return it.second;
}



VectorXd LineList::wave(std::string linestr) {
    int nn = nlines(linestr);

    if (nn==0) throw runtime_error("No such line available");
    VectorXd ret(nn);
    transform(llist.lower_bound(linestr), llist.upper_bound(linestr), &ret[0], _get_second);
    return ret;
}