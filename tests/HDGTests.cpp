
#define BOOST_TEST_MODULE hdg_unit_tests

#include <algorithm>
#include <complex>
#include <vector>
#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/math/constants/constants.hpp>


BOOST_AUTO_TEST_SUITE( linear_gs_operator_test_suite, * boost::unit_test::tolerance( 1e-14 ) )


BOOST_AUTO_TEST_CASE( _test )
{
	BOOST_TEST( foo == bar );
}
BOOST_AUTO_TEST_SUITE_END()
