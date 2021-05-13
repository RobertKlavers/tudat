/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_HYBRID_METHOD_MODEL_H
#define TUDAT_HYBRID_METHOD_MODEL_H

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/SimulationSetup/hybridOptimisationSettings.h"
#include <cmath>
#include <vector>
#include <Eigen/Dense>

namespace tudat
{
namespace low_thrust_trajectories
{

class HybridMethodModel
{
public:

    //! Constructor.
    HybridMethodModel( const Eigen::Vector6d& stateAtDeparture,
                     const Eigen::Vector6d& stateAtArrival,
                     const Eigen::VectorXd& initialCoStates,
                     const Eigen::VectorXd& finalCoStates,
                     const double maximumThrust,
                     const double specificImpulse,
                     const double timeOfFlight,
                     simulation_setup::NamedBodyMap& bodyMap,
                     const std::string bodyToPropagate,
                     const std::string centralBody,
                     std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
                     std::shared_ptr< simulation_setup::HybridOptimisationSettings > hybridOptimisationSettings ):
    stateAtDeparture_( stateAtDeparture ), stateAtArrival_( stateAtArrival ), initialCoStates_( initialCoStates ),
    finalCoStates_( finalCoStates ), maximumThrust_( maximumThrust ),
    specificImpulse_( specificImpulse ),
    timeOfFlight_( timeOfFlight ),
    bodyMap_( bodyMap ),
    bodyToPropagate_( bodyToPropagate ),
    centralBody_( centralBody ),
    integratorSettings_( integratorSettings ),
    hybridOptimisationSettings_( hybridOptimisationSettings )
    {
        // Initialise value of the total deltaV.
        totalDeltaV_ = 0.0;

        // Retrieve initial mass of the spacecraft.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Define function returning the current MEE costates.
        costatesFunction_ = [ = ]( const double currentTime )
        {
            Eigen::VectorXd currentCostates(initialCoStates.size());
            currentCostates.resize( 6 );

            for ( int i = 0 ; i < 6 ; i++ )
            {
                currentCostates[ i ] = initialCoStates_[ i ]
                        + ( currentTime / timeOfFlight_ ) * ( finalCoStates_[ i ] - initialCoStates_[ i ] );
            }
            return currentCostates;
        };

        // Initialise mass at time of flight (before propagation).
       Eigen::Vector6d propagatedStateAtTimeOfFlight = propagateTrajectory( );

    }


    //! Default destructor.
    ~HybridMethodModel( ) { }

    //! Retrieve MEE costates-based thrust acceleration.
    std::shared_ptr< simulation_setup::AccelerationSettings > getMEEcostatesBasedThrustAccelerationSettings( );

    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > getMEEcostatesBasedThrustMagnitudeSettings( );
    std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > getMEEcostatesBasedThrustDirectionSettings( );

    //! Retrieve hybrid method acceleration model (including thrust and central gravity acceleration)
    basic_astrodynamics::AccelerationMap getLowThrustTrajectoryAccelerationMap( );

    propagators::SingleArcDynamicsSimulator<> getDynamicsSimulator(
            double initialTime,
            double finalTime,
            Eigen::Vector6d initialState,
            double initialMass,
            std::shared_ptr<numerical_integrators::IntegratorSettings < double>> integratorSettings,
            bool withDependent = false);

    std::pair<std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd >> getTrajectoryOutput();

    //! Propagate the spacecraft trajectory to time of flight.
    Eigen::Vector6d propagateTrajectory( );

    //! Propagate the spacecraft trajectory to a given time.
    std::pair<Eigen::Vector6d, Eigen::Vector6d> propagateTrajectory( double initialTime, double finalTime, Eigen::Vector6d initialState, double initialMass);

    //! Propagate the trajectory to set of epochs as a function of theta.
    std::pair< double, Eigen::Vector6d > computeAverages(const Eigen::Vector6d& currentState, double currentTime, int numberOfSteps, double averagingTime );

    //! Propagate the trajectory to set of epochs.
    std::map< double, Eigen::Vector6d > propagateTrajectory(
            std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory );

    //! Return the deltaV associated with the thrust profile of the trajectory.
    double computeDeltaV( );

    //! Returns initial state at leg departure.
    Eigen::VectorXd getStateAtLegDeparture( )
    {
        return stateAtDeparture_;
    }

    //! Returns final state at leg arrival.
    Eigen::VectorXd getStateAtLegArrival( )
    {
        return stateAtArrival_;
    }

    //! Returns propagated mass when the time of flight is reached.
    double getMassAtTimeOfFlight( )
    {
        return massAtTimeOfFlight_;
    }

    //! Returns maximum allowed thrust.
    double getMaximumThrustValue( )
    {
        return maximumThrust_;
    }

    //! Returns time of flight.
    double getTimeOfFlight( )
    {
        return timeOfFlight_;
    }

    //! Return total deltaV required by the trajectory.
    double getTotalDeltaV( )
    {
        return computeDeltaV( );
    }

    //! Return the current MEE co-states.
    std::function< Eigen::VectorXd( const double ) > getCostatesFunction_( )
    {
        return costatesFunction_;
    }

    std::pair<std::vector<double>, Eigen::Vector6d> calculateFitness();

protected:

private:

    //! State vector of the vehicle at the leg departure.
    Eigen::Vector6d stateAtDeparture_;

    //! State vector of the vehicle at the leg arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Modified equinoctial elements.

    //! Initial co-states vector.
    Eigen::VectorXd initialCoStates_;

    //! Final co-states vector.
    Eigen::VectorXd finalCoStates_;

    //! Function returning the current MEE co-states.
    std::function< Eigen::VectorXd( const double ) > costatesFunction_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    double specificImpulse_;

    //! Time of flight for the leg.
    double timeOfFlight_;

    //! Pointer to the body map object.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    std::shared_ptr< simulation_setup::HybridOptimisationSettings > hybridOptimisationSettings_;

    //! Total deltaV.
    double totalDeltaV_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

    //! Mass of the spacecraft at the end of the propagation.
    double massAtTimeOfFlight_;

};


} // namespace low_thrust_trajectories
} // namespace tudat

#endif // TUDAT_HYBRID_METHOD_MODEL_H
