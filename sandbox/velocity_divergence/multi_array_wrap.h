#ifndef MULTI_ARRAY_WRAP_H
#define MULTI_ARRAY_WRAP_H


#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <algorithm>
#include <boost/static_assert.hpp>

namespace ma {

// Define some useful typedefs

typedef boost::multi_array_types::index_range r_;
typedef boost::multi_array_types::index i_;

// 4D arrays
// doubles
typedef boost::multi_array<double, 4> array4d;
typedef array4d::array_view<4>::type view4_4d;
typedef array4d::array_view<3>::type view4_3d;
typedef array4d::array_view<2>::type view4_2d;
typedef array4d::array_view<1>::type view4_1d;

// floats
typedef boost::multi_array<float, 4> array4f;
typedef array4f::array_view<4>::type  view4_4f;
typedef array4f::array_view<3>::type  view4_3f;
typedef array4f::array_view<2>::type  view4_2f;
typedef array4f::array_view<1>::type  view4_1f;

// 3D arrays
// doubles
typedef boost::multi_array<double, 3> array3d;
typedef array3d::array_view<3>::type  view3_3d;
typedef array3d::array_view<2>::type  view3_2d;
typedef array3d::array_view<1>::type  view3_1d;

// floats
typedef boost::multi_array<float, 3> array3f;
typedef array3f::array_view<3>::type  view3_3f;
typedef array3f::array_view<2>::type  view3_2f;
typedef array3f::array_view<1>::type  view3_1f;

// 2D arrays
// doubles
typedef boost::multi_array<double, 2> array2d;
typedef array2d::array_view<2>::type  view2_2d;
typedef array2d::array_view<1>::type  view2_1d;

// floats
typedef boost::multi_array<float, 2> array2f;
typedef array2f::array_view<2>::type  view2_2f;
typedef array2f::array_view<1>::type  view2_1f;



// Specialized loops
// Based on http://agentzlerich.blogspot.com/2010/01/providing-fill-and-foreach-algorithms.html

template<std::size_t dim>
struct _multi_for {
    BOOST_STATIC_ASSERT(dim > 1);

    // See http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for details about why MultiArray::iterator::reference is used below.
    template <class MultiArray, class Functor>
    void operator()(MultiArray& arr, Functor &func) {
        BOOST_STATIC_ASSERT(dim == MultiArray::dimensionality);
        _multi_for<dim-1> _f;
        for (typename MultiArray::iterator i = arr.begin(); i != arr.end(); ++i) {
            typename MultiArray::iterator::reference ri = *i;
            _f(ri, func);
        };
    }
};


template <>
struct _multi_for<1> {
    template <class MultiArray, class Functor>
    void operator()(MultiArray& arr, Functor &func) {
        for (typename MultiArray::iterator i = arr.begin(); i != arr.end(); ++i) {
            func(*i);
        };
    }
};

/* A simple nested loop iterator */
template <class MultiArray, class Functor>
void multi_for(MultiArray& arr, Functor& f) {
    _multi_for<MultiArray::dimensionality>() (arr, f);
}

////////////////////////////////////////////////////////////////////////////////////

template<std::size_t dim>
struct _multi_for_native_indices {
    BOOST_STATIC_ASSERT(dim > 1);

    // See http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for details about why MultiArray::iterator::reference is used below.
    template <class MultiArray, class Functor, class itype>
    void operator()(MultiArray& arr, Functor &func, itype& ilist, int idim) {
        BOOST_STATIC_ASSERT(dim == MultiArray::dimensionality);
        _multi_for_native_indices<dim-1> _f;
        int ii=0;
        for (typename MultiArray::iterator i = arr.begin() ; i != arr.end(); ++i, ++ii) {
            typename MultiArray::iterator::reference ri = *i;
            ilist[idim] = ii;
            _f(ri, func, ilist, idim+1);
        };
    }
};


template <>
struct _multi_for_native_indices<1> {
    template <class MultiArray, class Functor, class itype>
    void operator()(MultiArray& arr, Functor &func, itype& ilist, int idim) {
        int ii = 0;
        for (typename MultiArray::iterator i = arr.begin(); i != arr.end(); ++i, ++ii) {
            ilist[idim] = ii;
            func(*i, ilist);
        };
    }
};




/* A simple nested loop iterator */
template <class MultiArray, class Functor>
void multi_for_native_indices(MultiArray& arr, Functor& f) {
    boost::array<i_, MultiArray::dimensionality> ilist;
    _multi_for_native_indices<MultiArray::dimensionality>() (arr, f, ilist, 0);
}


////////////////////////////////////////////////////////////////////////////////////

template<std::size_t dim>
struct _multi_for_indices {
    BOOST_STATIC_ASSERT(dim > 1);

    // See http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for details about why MultiArray::iterator::reference is used below.
    template <class MultiArray, class Functor, class itype>
    void operator()(MultiArray& arr, Functor &func, itype& ilist, int idim) {
        BOOST_STATIC_ASSERT(dim == MultiArray::dimensionality);
        _multi_for_native_indices<dim-1> _f;
        int ii=arr.index_bases()[idim];
        for (typename MultiArray::iterator i = arr.begin() ; i != arr.end(); ++i, ++ii) {
            typename MultiArray::iterator::reference ri = *i;
            ilist[idim] = ii;
            _f(ri, func, ilist, idim+1);
        };
    }
};


template <>
struct _multi_for_indices<1> {
    template <class MultiArray, class Functor, class itype>
    void operator()(MultiArray& arr, Functor &func, itype& ilist, int idim) {
        int ii = arr.index_bases()[idim];
        for (typename MultiArray::iterator i = arr.begin(); i != arr.end(); ++i, ++ii) {
            ilist[idim] = ii;
            func(*i, ilist);
        };
    }
};




/* A simple nested loop iterator */
template <class MultiArray, class Functor>
void multi_for_indices(MultiArray& arr, Functor& f) {
    boost::array<i_, MultiArray::dimensionality> ilist;
    _multi_for_indices<MultiArray::dimensionality>() (arr, f, ilist, 0);
}





/**
 * Useful convenience functions.
 */


template<class T>
class _sum {
    private :
     T acc;
    public :
     _sum(T init) {acc = init;}
     void operator()(T& x) {acc += x;}
     T operator()() {return acc;}
};


template <class MultiArray, class T>
T sum(MultiArray& arr, T init) {
    _sum<T> f(init);
    multi_for(arr, f);
    return f();
}


} // Close Namespace ma
#endif // MULTI_ARRAY_WRAP_H
