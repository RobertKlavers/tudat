/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Propulsion/thrustGuidance.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace propulsion
{

//! Function to get the unit vector in the out-of-plane direction for increasing inclination
Eigen::Vector3d getForceDirectionOutOfPlane(
        const std::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const std::function< Eigen::Quaterniond() > rotationToGlobalFrameFunction,
        const double currentTime, const bool putForceInOppositeDirection )
{
    static Eigen::Vector6d currentState;
    currentStateFunction( currentState );

    Eigen::Matrix3d rotMat = rotationToGlobalFrameFunction().toRotationMatrix();

    Eigen::Vector3d bodyFixedThrustDirection = {0.0, 0.0, ( ( currentState[5] > 0 ) ? 1.0 : -1.0 )};

    Eigen::Vector3d inertialThrustVector = (rotMat * bodyFixedThrustDirection).normalized();

    Eigen::Vector3d currentVelocityDirection = currentState.segment( 3, 3).normalized( );

    std::cout << "\n----\n" << currentVelocityDirection << "\n----\n"  << rotMat << "\n----" << std::endl;

//
//     if v_z > 0 (currentstate[5]), thrust should be in the positive body-fixed up direction, down otherwise
//    Eigen::Vector6d outOfPlaneThrustDirection = (Eigen::Vector6d() << 0.0, 0.0, 0.0, 0.0, 0.0, ( ( currentState[5] > 0 ) ? 1.0 : -1.0 )).finished();
//    Eigen::Vector3d inertialThrustVector = transformStateToGlobalFrame( outOfPlaneThrustDirection, currentTime, rotationalEphemeris ).segment(3,3).normalized( );

//    std::cout << toInertialFrame << std::endl;
//
//    std::cout << "[" << currentTime << ", "
//              << currentState[0] << ", "
//              << currentState[1] << ", "
//              << currentState[2] << ", "
//              << inertialThrustVector[0] << ", "
//              << inertialThrustVector[1] << ", "
//              << inertialThrustVector[2] << "]" << std::endl;

    // Only use the velocity component of the inertial transformed body-fixed to inertial state
    return inertialThrustVector;
}



//! Function to get the unit vector colinear with velocity segment of a translational state.
Eigen::Vector3d getForceDirectionColinearWithVelocity(
        const std::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putForceInOppositeDirection )
{
    static Eigen::Vector6d currentState;
    currentStateFunction( currentState );

    Eigen::Vector3d currentPosition = currentState.segment(0,3);
    Eigen::Vector3d currentVelocity = currentState.segment(3,3);

    // Positive out-of-plane thrust for v_z > 0, negative for v_z < 0
    return currentPosition.cross(currentVelocity).normalized() * (currentVelocity[2] > 0 ? 1 : -1);


   // return ( ( putForceInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 3, 3 ) ).normalized( );
}

//! Function to get the unit vector colinear with position segment of a translational state.
Eigen::Vector3d getForceDirectionColinearWithPosition(
        const std::function< void( Eigen::Vector6d& ) > currentStateFunction,
        const double currentTime, const bool putForceInOppositeDirection )
{
    static Eigen::Vector6d currentState;
    currentStateFunction( currentState );
    return ( ( putForceInOppositeDirection == 1 ) ? -1.0 : 1.0 ) * ( currentState.segment( 0, 3 ) ).normalized( );
}

//! Function to get the force direction from a time-only function.
Eigen::Vector3d getForceDirectionFromTimeOnlyFunction(
        const double currentTime,
        const std::function< Eigen::Vector3d( const double ) > timeOnlyFunction )
{
    return timeOnlyFunction( currentTime ).normalized( );
}

} // namespace propulsion

} // namespace tudat

