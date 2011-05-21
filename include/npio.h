/* IO helper routines
 * Nikhil Padmanabhan, Yale, May, 2011
 *
 * The principal routine here is readAsciiFile, which can be specialized using different adaptors.
 * The current set of adaptors available are :
 *             - TupleAdaptor (to read into a boost::tuple
 *             - ArrayAdaptor (to read into *ANY* container that provides array (i.e. []) semantics).
 *
 *  Both of the routines are built around a simple LineProcessor
 */
#ifndef NPIO_H
#define NPIO_H

#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <sstream>

// Boost includes
// Tuples are nice things to read files with
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>

bool _isempty(const std::string &s) {return (s.size() == 0);}

std::list < std::string > tokenize(std::string ss, std::string delim) {
    std::list< std::string > strs;
    boost::split(strs, ss, boost::is_any_of(delim));
    strs.remove_if(_isempty);
    return strs;
}

template <class T, int N>
void LineProcessor(std::ifstream &ff, std::vector<T>& ll, std::string delim, bool (*parser)(T&, std::list< std::string >&)) {
    std::string ss;
    T elt;
    size_t pos;
    while (!ff.eof()) {
        getline(ff, ss);
        if (ss.empty()) continue;
        pos = ss.find_first_not_of(" \t");
        if ( (pos != std::string::npos) && (ss[0] != '#') ) {
            std::list< std::string > strs = tokenize(ss, delim) ;
            if (strs.size() != N) {continue;}
            if (parser(elt, strs)) ll.push_back(elt);
        }
    }
}

bool fill_tuple(const boost::tuples::null_type& , std::list < std::string >) {return true;}

template <class H, class T>
bool fill_tuple(boost::tuples::cons<H, T>& tup, std::list< std::string > &ll) {
    std::stringstream ss(ll.front());
    if ((ss >> tup.get_head()).fail()) {return false;}
    ll.pop_front();
    return fill_tuple(tup.get_tail(), ll);
}

template <class T>
bool tuple_parse(T& elt, std::list< std::string > &ll) {return fill_tuple(elt, ll);}

template <class T, int N>
void TupleAdaptor(std::ifstream &ff, std::vector<T>& ll, std::string delim) {
    LineProcessor<T, N>(ff, ll, delim, tuple_parse);
}

template <class T, int N>
bool array_parse(T& arr, std::list< std::string > &ll) {
    int ii = 0;
    BOOST_FOREACH (std::string s1, ll) {
        std::stringstream ss(s1);
        if ((ss >> arr[ii]).fail()) {return false;}
        ii++;
    }
    return true;
}


template <class T, int N>
void ArrayAdaptor(std::ifstream &ff, std::vector<T>& ll, std::string delim) {
    LineProcessor<T, N>(ff, ll, delim, array_parse<T,N>);
}

template <class T>
std::vector<T> readAsciiFile(std::string fn, void (*adaptor)(std::ifstream &, std::vector<T>&, std::string), std::string delim="\t ") {
    std::vector<T> ll;
    std::ifstream ff(fn.c_str());
    if (!ff.is_open()) {throw -999;} // Horrid design!
    adaptor(ff, ll, delim);
    ff.close();
    return ll;
}


#endif // NPIO_H
