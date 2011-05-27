#ifndef NP_FUNCTIONAL_H_
#define NP_FUNCTIONAL_H_

/* Some cool functional tricks */


#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>


template< typename R, typename I, typename S >
boost::function< R ( I ) > compose( boost::function< R ( I ) > f, boost::function< I ( S ) > g ) {
    return boost::lambda::bind( f, boost::lambda::bind( g, boost::lambda::_1 ) );
}

template< typename T, typename Y >
boost::function< T > operator *( boost::function< T > f, boost::function< Y > g ) {
    return compose( f, g );
}


#endif // NP_FUNCTIONAL_H_
