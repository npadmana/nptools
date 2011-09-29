#include <iostream>
#include <boost/phoenix.hpp>

#include <vector>
#include <algorithm>

using namespace std;
using namespace boost::phoenix;

template <typename F>
void print(F f) {
    cout << f() << endl;
}

int main() {
    using boost::phoenix::arg_names::arg1;


    cout << "Welcome to Phoenix \n";
    print(val(5));
    int i=10;
    print(ref(i));

    // Example
    // Allocate a vector
    vector<int> v(10, 0);
    // Now print it out
    for_each(v.begin(), v.end(), cout << arg1 << " ");
    cout << endl;

    // Now set to 1 to 10
    i = 0;
    for_each(v.begin(), v.end(), arg1=(ref(i)++));
    for_each(v.begin(), v.end(), cout << arg1 << " ");
    cout << endl;


    // Now do something relatively complicated
    // First, add 1
    // increase by 1 everything > 5
    for_each(v.begin(), v.end(), (
                 arg1++,
                 if_(arg1 > 5)[
                    arg1++
                 ]
                 ));
    for_each(v.begin(), v.end(), cout << arg1 << " ");
    cout << endl;

}
