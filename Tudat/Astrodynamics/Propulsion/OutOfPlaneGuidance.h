//
// Created by robert on 12-02-21.
//

#ifndef TUDATBUNDLE_OUTOFPLANEGUIDANCE_H
#define TUDATBUNDLE_OUTOFPLANEGUIDANCE_H

#include <cmath>

#include <functional>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/Propulsion/thrustGuidance.h"

namespace tudat
{

namespace propulsion
{

class OutOfPlaneGuidance: public BodyFixedForceDirectionGuidance
{
public:
    explicit OutOfPlaneGuidance(
            const std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
            const std::function< Eigen::Vector3d( ) > bodyFixedForceDirection =
            [ ]( ){ return  Eigen::Vector3d::UnitX( ); });

    ~OutOfPlaneGuidance();



private:
    //!  Function returning the state of the body under thrust as a function of time
    std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction_;

    void updateForceDirection( const double time );    //! Direction of thrust force in propagation frame, as computed by last call to updateForceDirection function.

    Eigen::Vector3d currentForceDirection_;

    Eigen::Vector3d getCurrentForceDirectionInPropagationFrame( )
    {
        return currentForceDirection_;
    }


    Eigen::Quaterniond getRotationToGlobalFrame( )
    {
        throw std::runtime_error( "Error, body-fixed frame to propagation frame not yet implemented for Out of Plane Guidance." );
    }

};

} // namespace propulsion

} // namespace tudat

#endif //TUDATBUNDLE_OUTOFPLANEGUIDANCE_H
