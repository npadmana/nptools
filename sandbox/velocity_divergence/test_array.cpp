#include "multi_array_wrap.h"
#include <iostream>
#include <boost/format.hpp>
#include <boost/static_assert.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

using namespace std;
using namespace ma;
using namespace boost;
using namespace boost::lambda;


template <class T> void index_test(double &x, T& ilist) {
    x = ilist[0] + ilist[1] + ilist[2];
}

template <class T> void index_test2(double &x, T& ilist) {
    x = ilist[0] + ilist[1] +10;
}


typedef array<i_, 3> itype;

void add_arr(double &x, itype ilist, array3d& b) {
    x += b(ilist);
}


int main() {

    cout << "Hello world\n";

    // Define a 3d array
    array3d a1(extents[3][3][3]);
    array3d b1(extents[3][3][3]);

    // Do a simple fill
    cout << "Filling array with 10s\n";
    multi_for(a1, _1 = 10.0);

    // Do a simple print
    multi_for(a1, cout << _1 << " ");
    cout << endl;

    // A more complicated fill
    double ii = 0.0;
    view3_2d a2 = a1[indices[1][r_()][r_()]];
    multi_for(a2, _1 = var(ii)++);

    // Do a simple print
    multi_for(a1, cout << _1 << " ");
    cout << endl;


    // Try a sum
    cout << "Sum of full array " << sum(a1, double(0.0)) << endl;
    cout << "Sum of subarray " << sum(a2, double(0.0)) << endl;


    // Do a more complicated fill
    multi_for_native_indices(a1, index_test<boost::array<i_, 3> >);
    multi_for(a1, cout << _1 << " ");
    cout << endl;

    multi_for_indices(a2, index_test2<boost::array<i_, 2> >);
    multi_for(a1, cout << _1 << " ");
    cout << endl;

    cout << "--------------\n";

    // Fill a1, and b1
    ii = 0.0;
    multi_for(a1, _1 = var(ii)++);
    multi_for(b1, _1 = var(ii)++);
    // Print these out
    multi_for(a1, cout << _1 << " ");
    cout << endl;
    multi_for(b1, cout << _1 << " ");
    cout << endl;

    // Now do something complicated with indices
    multi_for_native_indices(a1, bind(add_arr, _1, _2, var(b1)));
    // Print these out
    multi_for(a1, cout << _1 << " ");
    cout << endl;
    multi_for(b1, cout << _1 << " ");
    cout << endl;




}
