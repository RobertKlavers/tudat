//
// Created by robert on 12-02-21.
//

#include <cmath>

#include <functional>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Propulsion/OutOfPlaneGuidance.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"

namespace tudat {
namespace propulsion {

OutOfPlaneGuidance::OutOfPlaneGuidance(
        const std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
        const std::function< Eigen::Vector3d( ) > bodyFixedForceDirection)
        : BodyFixedForceDirectionGuidance( bodyFixedForceDirection ),
          thrustingBodyStateFunction_(thrustingBodyStateFunction) { }

void OutOfPlaneGuidance::updateForceDirection(const double time) {
    currentForceDirection_ = {0,0,1};
}

} // namespace propulsion

} // namespace tudat
