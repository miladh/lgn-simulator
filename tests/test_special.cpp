/**********************************************************************
 *  Test: spcial mathematical functions
 *
 *  Analytic source: by hand
 *
 * ********************************************************************/

#include <lgnSimulator.h>
#include <catch.hpp>


using namespace lgnSimulator;


TEST_CASE("second kind Bessel function", "[bessel]" ) {
    REQUIRE(Special::secondKindBessel(0.0) == 0);
    REQUIRE(Special::secondKindBessel(0.1) ==  Approx(0.049937526036241998).epsilon(1e-12));
    REQUIRE(Special::secondKindBessel(5.0) ==  Approx(-0.32757913759146522).epsilon(1e-12));
    REQUIRE(Special::secondKindBessel(7.01558666981562) ==  Approx(0.0).epsilon(1e-12));
    REQUIRE(Special::secondKindBessel(-7.01558666981562) ==  Approx(0.0).epsilon(1e-12));
    REQUIRE(Special::secondKindBessel(3.83170597020751) == Approx(0.0).epsilon(1e-12));
}


TEST_CASE("factorial function", "[factorial]" ) {
    REQUIRE(Special::factorial(0)==1);
    REQUIRE(Special::factorial(1)==1);
    REQUIRE(Special::factorial(4)==24);
    REQUIRE(Special::factorial(10)== 3628800);
    REQUIRE(Special::factorial(12)== 479001600);
}


TEST_CASE("heaviside function", "[heaviside]" ) {
    REQUIRE(Special::heaviside(-1.2)==0);
    REQUIRE(Special::heaviside(0.0)==1.);
    REQUIRE(Special::heaviside(2.2)==1.);
    REQUIRE(Special::heaviside(-100)==0);
    REQUIRE(Special::heaviside(1e-6)==1.);
    REQUIRE(Special::heaviside(34)==1.);
}


TEST_CASE("2d delta function", "[deltaVec2]" ) {
    REQUIRE(Special::delta({0,0},{0,0})==1);
    REQUIRE(Special::delta({-1e6, -1e6}, {-1e6, -1e6})==1);
    REQUIRE(Special::delta({3.456,3.456}, {3.456,3.456})==1);
    REQUIRE(Special::delta({2.1, 10.78},{2.1, 10.78})==1);
    REQUIRE(Special::delta({10.78, 2.1},{2.1, 10.78})==0);
    REQUIRE(Special::delta({-4.5, 900.89},{0.2, 67.6})==0);
    REQUIRE(Special::delta({2, -2},{2.0001, 2})==0);
}



TEST_CASE( "delta function", "[delta]" )
{
    REQUIRE(Special::delta(0,0) == 1);
    REQUIRE(Special::delta(-1e6,-1e6)== 1);
    REQUIRE(Special::delta(3.456,3.456)==1);
    REQUIRE(Special::delta(2.1,-2.1)==0);
    REQUIRE(Special::delta(1.,3.)==0);
    REQUIRE(Special::delta(-1230.,-50.)==0);

}


TEST_CASE("is odd function", "[isOddd]" ) {
    REQUIRE(Special::isOdd(0)==0);
    REQUIRE(Special::isOdd(-3)==1);
    REQUIRE(Special::isOdd(2)==0);
    REQUIRE(Special::isOdd(2.9)==0);
    REQUIRE(Special::isOdd(-3.6)==1);
    REQUIRE(Special::isOdd(91)==1);
}



TEST_CASE("nearest value", "[nearestValue]" ) {
    vec A = {0,1,2,3,4,5};
    double v_A = 2.3;
    REQUIRE(Special::nearestValue(A, v_A) == 2);

    vec B = {-2.3, 5.2, 8.2, 1.1, 2.67,-3.2};
    double v_B = -2.4;
    REQUIRE(Special::nearestValue(B, v_B) == -2.3);

    vec C = {9.66407446,  7.40204369,  8.47934683,  8.7268378 ,  9.537069  ,
             8.94007828,  6.37876932,  7.84503963,  8.70901142};
    double v_C = 7.5;
    REQUIRE(Special::nearestValue(C, v_C) == 7.40204369);


    vec D = { 4.20844744,  5.44088512, -1.44998235,  1.8764609 , -2.22633141,
              0.33623971,  7.23507673};
    double v_D = 0.0;
    REQUIRE(Special::nearestValue(D, v_D) == 0.33623971);
}





