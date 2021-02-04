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


#include <iostream>
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethodModel.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"

namespace tudat
{
namespace low_thrust_trajectories
{

using namespace orbital_element_conversions;

//! Retrieve MEE costates-based thrust acceleration.
std::shared_ptr< simulation_setup::AccelerationSettings > HybridMethodModel::getMEEcostatesBasedThrustAccelerationSettings( )
{
    // Define thrust direction settings from the MEE costates.
    std::shared_ptr< simulation_setup::MeeCostateBasedThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::MeeCostateBasedThrustDirectionSettings >( bodyToPropagate_, centralBody_,
                                                                                          costatesFunction_ );

    std::function< double ( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
      return specificImpulse_;
    };

    // Define bang-bang thrust magnitude settings based on MEE co-states.
    std::shared_ptr< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings >(
                maximumThrust_, specificImpulseFunction, costatesFunction_, bodyToPropagate_, centralBody_ );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    return thrustAccelerationSettings;

}


std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > HybridMethodModel::getMEEcostatesBasedThrustMagnitudeSettings( )
{
    std::function< double ( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
      return specificImpulse_;
    };

    // Return bang-bang thrust magnitude settings based on MEE co-states.
    return std::make_shared< simulation_setup::FromMeeCostatesBangBangThrustMagnitudeSettings >(
                maximumThrust_, specificImpulseFunction, costatesFunction_, bodyToPropagate_, centralBody_ );
}

std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > HybridMethodModel::getMEEcostatesBasedThrustDirectionSettings( )
{
    // Return thrust direction settings from the MEE costates.
    return std::make_shared< simulation_setup::MeeCostateBasedThrustDirectionSettings >( bodyToPropagate_, centralBody_, costatesFunction_ );
}


//! Retrieve hybrid method acceleration model (including thrust and central gravity acceleration)
basic_astrodynamics::AccelerationMap HybridMethodModel::getLowThrustTrajectoryAccelerationMap( )
{
    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
    bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getMEEcostatesBasedThrustAccelerationSettings( ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

    return accelerationModelMap;
}


//! Propagate the spacecraft trajectory to time-of-flight.
Eigen::Vector6d HybridMethodModel::propagateTrajectory( )
{
    Eigen::Vector6d propagatedState = propagateTrajectory( 0.0, timeOfFlight_, stateAtDeparture_, initialSpacecraftMass_ )[0];
    return propagatedState;
}



//! Propagate the spacecraft trajectory to a given time.
std::vector<Eigen::Vector6d>  HybridMethodModel::propagateTrajectory( double initialTime, double finalTime, Eigen::Vector6d initialState, double initialMass )
{
    // Re-initialise integrator settings.
    integratorSettings_->initialTime_ = initialTime;

    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ centralBody_ ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
//    bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getMEEcostatesBasedThrustAccelerationSettings( ) );

    // === CUSTOM TEST THRUST!  ===
//    /*
    double thrustMagnitude = 10.0;
    double specificImpulse = 2000.0;
    std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< simulation_setup::ThrustDirectionFromStateGuidanceSettings >( centralBody_, true, false );
    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< simulation_setup::ConstantThrustMagnitudeSettings >( thrustMagnitude, specificImpulse );

    // Define acceleration model settings.
    bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back(
            std::make_shared< simulation_setup::ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings ) );
//    */

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );


    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationModelMap );

    // Ensure that the propagation stops when the required time of flight is required.
    std::shared_ptr< propagators::PropagationTimeTerminationSettings > terminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( finalTime, true );


        // Define propagator settings.
        std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > translationalStatePropagatorSettings =
                std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >{ centralBody_ }, accelerationModelMap, std::vector< std::string >{ bodyToPropagate_ },
                    initialState, terminationSettings, propagators::gauss_modified_equinoctial );

        // Create settings for propagating the mass of the vehicle.
        std::shared_ptr< propagators::MassPropagatorSettings< double > > massPropagatorSettings
                = std::make_shared< propagators::MassPropagatorSettings< double > >(
                    std::vector< std::string >{ bodyToPropagate_ }, massRateModel,
                    ( Eigen::Matrix< double, 1, 1 >( ) << bodyMap_[ bodyToPropagate_ ]->getBodyMass( ) ).finished( ),
                    terminationSettings );

        // Create list of propagation settings.
        std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
        propagatorSettingsVector.push_back( translationalStatePropagatorSettings );
        propagatorSettingsVector.push_back( massPropagatorSettings );

        // Hybrid propagation settings.
        std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings =
                std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

        integratorSettings_->initialTime_ = initialTime;

        // Propagate the trajectory.
        propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap_, integratorSettings_, propagatorSettings );

        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );

        double propagationResultTime = numericalSolution.rbegin( )->first;
        Eigen::VectorXd propagationResult = numericalSolution.at(propagationResultTime);

        //TODO SEE WHAT HAPPENS
//        std::cout << " -- computeStateDerivative t: " << propagationResultTime << " -> \n " << propagationResult << std::endl;

        double gravitationalParameter = bodyMap_[centralBody_]->getGravityFieldModel()->getGravitationalParameter();

        Eigen::Vector6d propagatedMEEState = propagationResult.segment(0, 6);
        Eigen::Vector6d computedCartesianState = orbital_element_conversions::convertModifiedEquinoctialToCartesianElements(propagatedMEEState, gravitationalParameter, false);


        // Find the state derivatives at t_f
        Eigen::VectorXd currentStateDerivative;
        currentStateDerivative = dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(
                propagationResultTime, propagationResult);

        Eigen::Vector6d computedMMEStateDerivatives = currentStateDerivative.segment(0, 6);

        std::cout << "-- computedMMEStateDerivatives -- \n" << computedMMEStateDerivatives << std::endl;
        std::cout << "-- computedCartesianState -- \n" << computedCartesianState << std::endl;

    // Retrieve state and mass of the spacecraft at the end of the propagation.
    if ( finalTime == timeOfFlight_ )
    {
        massAtTimeOfFlight_ = propagationResult[ 6 ];
    }

    return std::vector<Eigen::Vector6d> {computedCartesianState, computedMMEStateDerivatives};
}


    Eigen::Vector6d HybridMethodModel::propagateTrajectoryForTheta(
            std::map< double, Eigen::Vector6d >& propagatedTrajectory, int numberOfSteps )
    {
        // Initialise propagated state.
        Eigen::Vector6d propagatedState = stateAtDeparture_;

        std::vector<Eigen::Vector6d> propagationResult;
        std::vector<Eigen::Vector6d> stateDerivatives;

        // Generate Epochs for Theta
        std::vector< double > thetaEpochs;
        double deltaTheta = (2.0 * mathematical_constants::PI) / numberOfSteps;

        for ( int i = 1 ; i <= numberOfSteps ; i++ )
        {
            thetaEpochs.push_back( deltaTheta * i );
        }

        // Initialise mass of the spacecraft at departure.
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
        double currentMass = initialSpacecraftMass_;

        double gravitationalParameter = bodyMap_[centralBody_]->getGravityFieldModel()->getGravitationalParameter();

        // Temp Debug print
        Eigen::Vector6d keplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(propagatedState,
                                                                                                         gravitationalParameter);
        std::cout << "-- Keplerian Elements -- \n " << keplerianState << std::endl;

        // Convert from theta to time (container vector)
        std::vector< double > timeEpochs;

        double currentTime = 0.0;

        for ( int epochIndex = 0 ; epochIndex < thetaEpochs.size( ) ; epochIndex++ )
        {
            double currentTheta = thetaEpochs[ epochIndex ];
            timeEpochs.push_back(currentTime);

            double a = orbital_element_conversions::convertCartesianToKeplerianElements(propagatedState,
                                                                             gravitationalParameter)[0];
            double r = propagatedState.segment(0 ,3 ).norm();

            // DEBUG
            std::cout << "r: " << r << ", a: " << a << std::endl;

            // Convert to time delta
            double deltaTime = deltaTheta *
                               (r / (a * std::sqrt(gravitationalParameter / std::pow(a, 3.0))));


            std::cout << "Epoch: " << epochIndex << ", cTime: " << currentTime << ", cTheta: " << currentTheta << ", dT: " << deltaTime << std::endl;

            if ( epochIndex == 0 )
            {
                if ( currentTheta > 0.0 )
                {
                    propagationResult = propagateTrajectory( currentTime, currentTime + deltaTime, propagatedState, currentMass );

                    currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
                    propagatedState = propagationResult[0];

                    stateDerivatives.push_back(propagationResult[1]);
                }
            }
            else
            {
                propagationResult = propagateTrajectory( currentTime, currentTime + deltaTime, propagatedState, currentMass );

                currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
                propagatedState = propagationResult[0];

                stateDerivatives.push_back(propagationResult[1]);
            }
            propagatedTrajectory[ currentTime ] =  propagatedState;
            currentTime = currentTime + deltaTime;
        }

        double orbitPeriod = currentTime;

        bodyMap_[ centralBody_ ]->setConstantBodyMass( initialSpacecraftMass_ );

        // Create Trapezoidal Quadrature Integrator
        tudat::numerical_quadrature::TrapezoidNumericalQuadrature< double, Eigen::Vector6d > integrator(
                timeEpochs, stateDerivatives );
        Eigen::Vector6d computedIntegralTrapezoid = integrator.getQuadrature();
        std::cout << "-- Trapezoidal Integral: \n" << computedIntegralTrapezoid << std::endl;

        Eigen::Vector6d averageStateProgression = computedIntegralTrapezoid / orbitPeriod;

        return averageStateProgression;
    }
//! Propagate the trajectory to set of epochs.
std::map< double, Eigen::Vector6d > HybridMethodModel::propagateTrajectory(
        std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory )
{
    // Initialise propagated state.
    Eigen::Vector6d propagatedState = stateAtDeparture_;

    // Initialise mass of the spacecraft at departure.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );
    double currentMass = initialSpacecraftMass_;


    for ( int epochIndex = 0 ; epochIndex < epochs.size( ) ; epochIndex++ )
    {
        double currentTime = epochs[ epochIndex ];
        if ( epochIndex > 0 )
        {
            if ( currentTime < epochs[ epochIndex - 1 ] )
            {
                throw std::runtime_error( "Error when propagating trajectory with hybrid method, epochs at which the trajectory should be "
                                          "computed are not in increasing order." );
            }
        }
        if ( ( currentTime < 0.0 ) || ( currentTime > timeOfFlight_ ) )
        {
            throw std::runtime_error( "Error when propagating trajectory with hybrid method, epochs at which the trajectory should be "
                                      "computed are not constrained between 0.0 and timeOfFlight." );
        }


        if ( epochIndex == 0 )
        {
            if ( currentTime > 0.0 )
            {
                propagatedState = propagateTrajectory( 0.0, currentTime, propagatedState, currentMass )[0];
                currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
            }
            propagatedTrajectory[ currentTime ] = propagatedState;
        }
        else
        {
            propagatedState = propagateTrajectory( epochs[ epochIndex - 1 ], currentTime, propagatedState, currentMass )[0];
            currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
            propagatedTrajectory[ currentTime ] = propagatedState;
        }

    }

    bodyMap_[ centralBody_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    return propagatedTrajectory;
}


//! Return the deltaV associated with the thrust profile of the trajectory.
double HybridMethodModel::computeDeltaV( )
{

    // Compute (constant) mass rate.
    double massRate = - maximumThrust_ /
            ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );

    // Compute time during which the engine was switched on.
    double engineSwitchedOnDuration = ( massAtTimeOfFlight_ - initialSpacecraftMass_ ) / massRate;

    std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_ =
            std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( 0.0, 16 );

    // Thrust acceleration function to use quadrature.
    // Define thrust acceleration as a function of time (to be integrated to compute the associated deltaV).
    std::function< double( const double ) > thrustAcceleration = [ = ] ( const double currentTime ){

        // Compute current mass.
        double currentMass = initialSpacecraftMass_ + massRate * currentTime;

        // Compute and return current thrust acceleration.
        double currentThrustAcceleration = maximumThrust_ / currentMass;

        return currentThrustAcceleration;
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( thrustAcceleration, quadratureSettings_, engineSwitchedOnDuration );

    // Compute deltaV analytically.
    double deltaV = - specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION
            * std::log( 1.0 + massRate / initialSpacecraftMass_ * engineSwitchedOnDuration );

    return deltaV;
}



} // namespace low_thrust_trajectories
} // namespace tudat
