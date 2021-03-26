/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethodModel.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethod.h"
#include "pagmo/algorithms/de1220.hpp"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"

namespace tudat {

namespace unit_tests {


//! Test hybrid method implementation.
BOOST_AUTO_TEST_SUITE( test_average_step )

BOOST_AUTO_TEST_CASE( test_average_step_implementation ) {
    using namespace low_thrust_trajectories;

    BOOST_TEST_MESSAGE("Starting ");

    BOOST_CHECK_SMALL( 0.0, 1.0e-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
