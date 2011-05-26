#ifndef MULTI_ARRAY_WRAP_H
#define MULTI_ARRAY_WRAP_H


#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace ma {

using namespace boost::lambda;

// Define some useful typedefs

typedef boost::multi_array_types::index_range r_;
typedef boost::multi_array_types::index i_;


template <class Array>
class MA {
   public :
    typedef typename Array::value_type value_type;
    typedef typename Array::reference reference;
    typedef typename Array::size_type size_type;
    typedef typename Array::difference_type difference_type;
    typedef typename Array::iterator iterator;
    typedef typename Array::index index;
    typedef typename Array::element element;
    static const size_type dimensionality = Array::dimensionality;

   private :
    boost::shared_ptr<Array> ptr;

   public :
    // Simple constructor
     MA(Array * p) : ptr(p) {}


     // Basic Accessor functions
     const size_type* shape() {return ptr->shape();}
     const index* strides() {return ptr->strides();}
     const index* index_bases() {return ptr->index_bases();}
     const element* origin() {return ptr->a.origin();}
     size_type num_dimensions() {return ptr->num_dimensions();}
     size_type size() {return ptr->size();}
     template <class IndexList> element& operator()(const IndexList& list ) {return (*ptr)(list);}
     template <class IndexList> element& operator()(const IndexList& list ) const {return (*ptr)(list);}
     iterator begin() {return ptr->begin();}
     iterator end() {return ptr->end();}
     Array& operator()() {return *ptr;}
     reference operator[](index i) { return (*ptr)[i];}
     MA<reference> sub(index i) {
         return MA<reference> ( new reference((*ptr)[i]) ) ;
     }
     template <class A, class RangeList> MA<A> view(const RangeList& list) {
         return MA<A> ( new A((*ptr)[list]) ) ;
     }
     element* ref() {
         boost::array<i_, dimensionality> tmp;
         // Correct for strides
         std::transform(index_bases(), index_bases() + dimensionality,
                        strides(), tmp.begin(), std::multiplies<index>());
         // add the sum of the products to the origin
         return std::accumulate(tmp.begin(), tmp.end(), origin());
     }

     // Simple arithmetic overloads
     void set(element x) { multi_for(*ptr, _1=x);}

     // Note that we use the fact that we have a lightweight copy of MA objects
     void operator+=(element x) { multi_for(*ptr, _1 += x);}
     template <class A> void operator+=(MA<A> x) { multi_for(*ptr, x, _1 += _2);}
     void operator-=(element x) { multi_for(*ptr, _1 -= x);}
     template <class A> void operator-=(MA<A> x) { multi_for(*ptr, x, _1 -= _2);}
     void operator*=(element x) { multi_for(*ptr, _1 *= x);}
     template <class A> void operator*=(MA<A> x) { multi_for(*ptr, x, _1 *= _2);}
     void operator/=(element x) { multi_for(*ptr, _1 /= x);}
     template <class A> void operator/=(MA<A> x) { multi_for(*ptr, x, _1 /= _2);}
     void operator%=(element x) { multi_for(*ptr, _1 %= x);}
     template <class A> void operator%=(MA<A> x) { multi_for(*ptr, x, _1 %= _2);}



};


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


// 1D arrays
typedef boost::multi_array<double, 1> array1d;
typedef array1d::array_view<1>::type view1_1d;
typedef boost::multi_array<float, 1> array1f;
typedef array1f::array_view<1>::type view1_1f;
typedef boost::multi_array<int, 1> array1i;
typedef array1i::array_view<1>::type view1_1i;

/////////////////////////////////////////////////////////////////////////////////
template <class A, class B>
bool compatible(A& x, B& y) {
    if (A::dimensionality != B::dimensionality) return false;
    for (int ii=0; ii<A::dimensionality; ++ii)
        if (x.shape()[ii] > y.shape()[ii]) return false;
    return true;
}




////////////////////////////////////////////////////////////////////////////////////

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

// Specialized loops
// Based on http://agentzlerich.blogspot.com/2010/01/providing-fill-and-foreach-algorithms.html

template<std::size_t dim>
struct _multi_for2 {
    BOOST_STATIC_ASSERT(dim > 1);

    // See http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for details about why MultiArray::iterator::reference is used below.
    template <class MA1, class MA2, class Functor>
    void operator()(MA1& arr1, MA2 &arr2, Functor &func) {
        BOOST_STATIC_ASSERT(dim == MA1::dimensionality);
        BOOST_STATIC_ASSERT(dim == MA2::dimensionality);
        _multi_for2<dim-1> _f;
        typename MA2::iterator i2 = arr2.begin();
        for (typename MA1::iterator i1 = arr1.begin(); i1 != arr1.end(); ++i1, ++i2) {
            typename MA1::iterator::reference ri1 = *i1;
            typename MA2::iterator::reference ri2 = *i2;
            _f(ri1, ri2, func);
        };
    }
};


template <>
struct _multi_for2<1> {
    template <class MA1, class MA2, class Functor>
    void operator()(MA1& arr1, MA2& arr2, Functor &func) {
        typename MA2::iterator i2 = arr2.begin();
        for (typename MA1::iterator i1 = arr1.begin(); i1 != arr1.end(); ++i1, ++i2) {
            func(*i1, *i2);
        };
    }
};

/* A simple nested loop iterator */
template <class MA1, class MA2, class Functor>
void multi_for(MA1& arr1, MA2& arr2, Functor& f) {
    BOOST_STATIC_ASSERT(MA1::dimensionality == MA2::dimensionality);
    if (!compatible(arr1, arr2) ) throw std::out_of_range("Arrays 1 and 2 are not compatible");
    _multi_for2<MA1::dimensionality>() (arr1, arr2, f);
}

////////////////////////////////////////////////////////////////////////////////////

// Specialized loops
// Based on http://agentzlerich.blogspot.com/2010/01/providing-fill-and-foreach-algorithms.html

template<std::size_t dim>
struct _multi_for3 {
    BOOST_STATIC_ASSERT(dim > 1);

    // See http://groups.google.com/group/boost-list/browse_thread/thread/e16f32c4411dea08
    // for details about why MultiArray::iterator::reference is used below.
    template <class MA1, class MA2, class MA3, class Functor>
    void operator()(MA1& arr1, MA2 &arr2, MA3 &arr3, Functor &func) {
        BOOST_STATIC_ASSERT(dim == MA1::dimensionality);
        BOOST_STATIC_ASSERT(dim == MA2::dimensionality);
        BOOST_STATIC_ASSERT(dim == MA3::dimensionality);
        _multi_for3<dim-1> _f;
        typename MA2::iterator i2 = arr2.begin();
        typename MA3::iterator i3 = arr3.begin();
        for (typename MA1::iterator i1 = arr1.begin(); i1 != arr1.end(); ++i1, ++i2, ++i3) {
            typename MA1::iterator::reference ri1 = *i1;
            typename MA2::iterator::reference ri2 = *i2;
            typename MA3::iterator::reference ri3 = *i3;
            _f(ri1, ri2, ri3, func);
        };
    }
};


template <>
struct _multi_for3<1> {
    template <class MA1, class MA2, class MA3, class Functor>
    void operator()(MA1& arr1, MA2& arr2, MA3& arr3, Functor &func) {
        typename MA2::iterator i2 = arr2.begin();
        typename MA3::iterator i3 = arr3.begin();
        for (typename MA1::iterator i1 = arr1.begin(); i1 != arr1.end(); ++i1, ++i2, ++i3) {
            func(*i1, *i2, *i3);
        };
    }
};

/* A simple nested loop iterator */
template <class MA1, class MA2, class MA3, class Functor>
void multi_for(MA1& arr1, MA2& arr2, MA3& arr3, Functor& f) {
    BOOST_STATIC_ASSERT(MA1::dimensionality == MA2::dimensionality);
    BOOST_STATIC_ASSERT(MA1::dimensionality == MA3::dimensionality);
    if (!compatible(arr1, arr2) ) throw std::out_of_range("Arrays 1 and 2 are not compatible");
    if (!compatible(arr1, arr3) ) throw std::out_of_range("Arrays 1 and 3 are not compatible");
    _multi_for3<MA1::dimensionality>() (arr1, arr2, arr3, f);
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




/***
  * Operator helper functions
  */





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
