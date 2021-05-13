/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *
 *    Kluever (2010), Low-Thrust Trajectory Optimization Using Orbital Averaging and Control Parameterization, In: Conway,
 *    (editor) Spacecraft trajectory optimization. Cambridge University Press, 2010.
 *    Boudestijn (2014), DEVELOPMENT OF A LOW -THRUST EARTH-CENTERED TRANSFER OPTIMIZER FOR THE PRELIMINARY MISSION DESIGN PHASE,
 *    M.Sc. Thesis, Delft University of Technology
 */

#ifndef TUDAT_HYBRID_OPTIMISATION_SETTINGS_H
#define TUDAT_HYBRID_OPTIMISATION_SETTINGS_H

#include <Eigen/Geometry>
#include <boost/bind.hpp>
#include <functional>
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"
#include "pagmo/algorithm.hpp"


namespace tudat {

    namespace simulation_setup {
        class HybridOptimisationSettings {
        public:
            HybridOptimisationSettings(const Eigen::Vector6d epsilonUpper, const Eigen::Vector6d constraintWeights,
                                       const double weightMass, const double weightTimeOfFlight, const bool debug) :
                    epsilonUpper_(epsilonUpper),
                    constraintWeights_(constraintWeights),
                    weightMass_(weightMass),
                    weightTimeOfFlight_(weightTimeOfFlight), debug_(debug) {}

            //! Destructor.
            virtual ~HybridOptimisationSettings() {}

            const Eigen::Vector6d epsilonUpper_;
            const Eigen::Vector6d constraintWeights_;
            const double weightMass_;
            const double weightTimeOfFlight_;
            const bool debug_;
        };


    } // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_HYBRID_OPTIMISATION_SETTINGS_H
