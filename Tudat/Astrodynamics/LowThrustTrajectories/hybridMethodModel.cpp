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

propagators::SingleArcDynamicsSimulator<> HybridMethodModel::getDynamicsSimulator(
        double initialTime,
        double finalTime,
        Eigen::Vector6d initialState,
        double initialMass,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        bool withDependent
        ) {
    // Re-initialise integrator settings.
    integratorSettings->initialTime_ = initialTime;

    bodyMap_[bodyToPropagate_]->setConstantBodyMass(initialMass);

    // std::cout << "t0: " << initialTime << ", tf: " << finalTime << std::endl;

    // // === CUSTOM TEST THRUST!  ===
    // TODO:REMOVE OR REPLACE
    // std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > colinearThrustDirectionGuidanceSettings =
    //         std::make_shared< simulation_setup::ThrustDirectionFromStateGuidanceSettings >( centralBody_, true, false );
    // std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > constantThrustMagnitudeSettings =
    //         std::make_shared< simulation_setup::ConstantThrustMagnitudeSettings >( maximumThrust_, specificImpulse_ );
    // std::shared_ptr< simulation_setup::AccelerationSettings > constantThrustAccelerationSettings = std::make_shared< simulation_setup::ThrustAccelerationSettings >( colinearThrustDirectionGuidanceSettings, constantThrustMagnitudeSettings );

    basic_astrodynamics::AccelerationMap accelerationModelMap = getLowThrustTrajectoryAccelerationMap();

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap_, accelerationModelMap );

    // Ensure that the propagation stops when the required time of flight is required.
    std::shared_ptr< propagators::PropagationTimeTerminationSettings > terminationSettings
            = std::make_shared< propagators::PropagationTimeTerminationSettings >( finalTime, true );

    std::vector<std::shared_ptr<propagators::SingleDependentVariableSaveSettings> > dependentVariablesList;

    dependentVariablesList.push_back(std::make_shared<propagators::SingleDependentVariableSaveSettings>(
            propagators::keplerian_state_dependent_variable, bodyToPropagate_, centralBody_));
    dependentVariablesList.push_back(std::make_shared<propagators::SingleAccelerationDependentVariableSaveSettings>(
            basic_astrodynamics::thrust_acceleration, bodyToPropagate_, bodyToPropagate_, 0));
    dependentVariablesList.push_back(std::make_shared<propagators::SingleDependentVariableSaveSettings>(
            propagators::lvlh_to_inertial_frame_rotation_dependent_variable, bodyToPropagate_, centralBody_));

    // dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
    //         propagators::rotation_matrix_to_body_fixed_frame_variable, "Vehicle", "Earth" ) );
    std::shared_ptr<propagators::DependentVariableSaveSettings> dependentVariablesToSave =
            std::make_shared<propagators::DependentVariableSaveSettings>(dependentVariablesList, withDependent);

    // Define propagator settings.
    std::shared_ptr<propagators::TranslationalStatePropagatorSettings<double> > translationalStatePropagatorSettings =
            std::make_shared<propagators::TranslationalStatePropagatorSettings<double> >(
                    std::vector<std::string>{centralBody_}, accelerationModelMap,
                    std::vector<std::string>{bodyToPropagate_},
                    initialState, terminationSettings, propagators::gauss_modified_equinoctial);

    // Create settings for propagating the mass of the vehicle.
    std::shared_ptr<propagators::MassPropagatorSettings<double> > massPropagatorSettings
            = std::make_shared<propagators::MassPropagatorSettings<double> >(
                    std::vector<std::string>{bodyToPropagate_}, massRateModel,
                    (Eigen::Matrix<double, 1, 1>() << bodyMap_[bodyToPropagate_]->getBodyMass()).finished(),
                    terminationSettings);

    // Create list of propagation settings.
    std::vector<std::shared_ptr<propagators::SingleArcPropagatorSettings<double> > > propagatorSettingsVector;
    propagatorSettingsVector.push_back(translationalStatePropagatorSettings);
    propagatorSettingsVector.push_back(massPropagatorSettings);

    // Hybrid propagation settings.
    std::shared_ptr<propagators::PropagatorSettings<double> > propagatorSettings =
            std::make_shared<propagators::MultiTypePropagatorSettings<double> >(propagatorSettingsVector,
                                                                                terminationSettings,
                                                                                dependentVariablesToSave);

    // Propagate the trajectory.
    propagators::SingleArcDynamicsSimulator<> dynamicsSimulator(bodyMap_, integratorSettings, propagatorSettings);

    return dynamicsSimulator;
}


//! Propagate the spacecraft trajectory to time-of-flight.
Eigen::Vector6d HybridMethodModel::propagateTrajectory( )
{
    if (hybridOptimisationSettings_->debug_) {
        std::cout << "    |" << timeOfFlight_ << "," << stateAtDeparture_.transpose()<< "," << initialSpacecraftMass_ << std::endl;
    }
    Eigen::Vector6d propagatedState = propagateTrajectory( 0.0, timeOfFlight_, stateAtDeparture_, initialSpacecraftMass_ ).first;
    return propagatedState;
}

//! Propagate the spacecraft trajectory to a given time.
std::pair<Eigen::Vector6d, Eigen::Vector6d> HybridMethodModel::propagateTrajectory( double initialTime, double finalTime, Eigen::Vector6d initialState, double initialMass)
{
    integratorSettings_->initialTime_ = initialTime;

    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator = getDynamicsSimulator(initialTime, finalTime, initialState, initialMass, integratorSettings_);

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
    Eigen::VectorXd propagationOriginalResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second;

    // Retrieve state and mass of the spacecraft at the end of the propagation.
    Eigen::Vector6d propagatedState = propagationOriginalResult.segment( 0, 6 );
    Eigen::Vector6d computedMMEStateDerivatives;

    if ( finalTime == timeOfFlight_ )
    {
        massAtTimeOfFlight_ = propagationOriginalResult[ 6 ];
    }

    return std::pair<Eigen::Vector6d, Eigen::Vector6d> (propagatedState, computedMMEStateDerivatives);
}

//! Propagate the trajectory to set of epochs.
std::map< double, Eigen::Vector6d > HybridMethodModel::propagateTrajectory(
        std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory )
{
    std::cout<< "PROPAGATING TO EPOCHS SHOULD BE OBSOLETE? " << std::endl;
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
                propagatedState = propagateTrajectory( 0.0, currentTime, propagatedState, currentMass).first;
                currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
            }
            propagatedTrajectory[ currentTime ] = propagatedState;
        }
        else
        {
            propagatedState = propagateTrajectory( epochs[ epochIndex - 1 ], currentTime, propagatedState, currentMass).first;
            currentMass = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
            propagatedTrajectory[ currentTime ] = propagatedState;
        }

    }

    bodyMap_[ centralBody_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    return propagatedTrajectory;
}

//! Utility to retrieve the integration result and all saved dependent variables
std::pair<std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd >> HybridMethodModel::getTrajectoryOutput() {
    std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings =
            std::make_shared<numerical_integrators::IntegratorSettings<double> >
                    (numerical_integrators::rungeKutta4, 0.0, 60.0);

    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator = getDynamicsSimulator(0.0, timeOfFlight_, stateAtDeparture_, initialSpacecraftMass_, integratorSettings, true);
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    // index to denote where the thrust accelerations start in the dependent variable history.
    int thrustIndex = 6;
    // Manually add thrust force in LVLH frame to output
    for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
         outputIterator != dependentVariableResult.end( ); outputIterator++ )
    {
        Eigen::Matrix3d currentRotationMatrix =
                propagators::getMatrixFromVectorRotationRepresentation( outputIterator->second.segment( thrustIndex+3, 9 ) );
        Eigen::Vector3d currentThrust = outputIterator->second.segment( thrustIndex+0, 3 );
        Eigen::VectorXd newOutput = Eigen::VectorXd( thrustIndex+15 );
        newOutput.segment( 0, thrustIndex+12 ) = outputIterator->second;
        newOutput.segment( thrustIndex+12, 3 ) =
                integrationResult.at( outputIterator->first )( 6 ) *
                ( currentRotationMatrix.transpose( ) * currentThrust );
        dependentVariableResult[ outputIterator->first ] = newOutput;
    }

    return std::pair<std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd >> (integrationResult, dependentVariableResult);
}

std::pair<Eigen::VectorXd, Eigen::Vector6d> HybridMethodModel::calculateFitness() {
    // Propagate until time of flight is reached.
    Eigen::Vector6d finalPropagatedState = propagateTrajectory( );

    // Convert final propagated state to MEE.
    Eigen::Vector6d finalPropagatedKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
            finalPropagatedState, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );

    // Convert targeted final state to MEE.
    Eigen::Vector6d finalTargetedKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
            stateAtArrival_, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );


    // TODO: Wrap Angle? [0, 360)
    Eigen::Vector6d error = (finalPropagatedKeplerianState - finalTargetedKeplerianState).cwiseAbs();

    // Compute auxiliary variables.
    Eigen::Vector6d c;
    Eigen::Vector6d r;
    Eigen::Vector6d epsilon;
    Eigen::Vector6d epsilon_lower = Eigen::Vector6d::Zero(6);
    Eigen::Vector6d epsilon_upper = Eigen::Vector6d::Zero(6);

    // TODO this is ugly, make this more elegant
    double deg2rad = mathematical_constants::PI / 180.0;
    epsilon_upper(0) = hybridOptimisationSettings_->epsilonUpper_(0),
    epsilon_upper(1) = hybridOptimisationSettings_->epsilonUpper_(1),
    epsilon_upper(2) = hybridOptimisationSettings_->epsilonUpper_(2) * deg2rad,
    epsilon_upper(3) = hybridOptimisationSettings_->epsilonUpper_(3) * deg2rad,
    epsilon_upper(4) = hybridOptimisationSettings_->epsilonUpper_(4) * deg2rad,
    epsilon_upper(5) = hybridOptimisationSettings_->epsilonUpper_(5) * deg2rad;

    // Aggregate Objective Function (6 epsilons + tof + mass)
    Eigen::VectorXd aof_vector(8);

    // 5 for no longitude targeting?
    for ( int i = 0 ; i < 6 ; i++ )
    {
        c[ i ] = 1.0 / (epsilon_upper(i)- epsilon_lower(i));
        r[ i ] = 1.0 - epsilon_upper(i) * c[ i ];
        epsilon[ i ] = error[ i ] * c[ i ] + r[ i ];
    }
    // Determine scaled epsilon at TOF
    for ( int i = 0 ; i < 6 ; i++) {
        aof_vector(i) = (hybridOptimisationSettings_->constraintWeights_[i] * epsilon[i] * epsilon[i]);
    }

    // AOF From Jimenez: F = sum(orbit error) + W_t*t_f + W_m*(1 - m_f/m_0)
    double finalMass = getMassAtTimeOfFlight();
    aof_vector(6) = (hybridOptimisationSettings_->weightTimeOfFlight_ * timeOfFlight_/physical_constants::JULIAN_DAY);
    aof_vector(7) = (hybridOptimisationSettings_->weightMass_*(1 - finalMass/initialSpacecraftMass_));

    // Some debugging statements we don't actually directly output the fitness
    if (hybridOptimisationSettings_->debug_) {
        double totalEps = aof_vector.sum();
        std::cout << "\n\n--fitcalc--" << std::endl;
        std::cout << "  cfu:" << costatesFunction_(0.0).transpose() << "] | [" << costatesFunction_(timeOfFlight_).transpose() << std::endl;
        std::cout << "  x_f: " << finalPropagatedState.transpose() << std::endl;
        std::cout << "  err:" << error.transpose() << std::endl;
        std::cout << "  sum_e: " << totalEps << std::endl;
    }
    return {aof_vector, error};
}

//! Return the deltaV associated with the thrust profile of the trajectory.
double HybridMethodModel::computeDeltaV( )
{
    // Compute deltaV analytically.
    return specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION
                     * std::log( initialSpacecraftMass_ / massAtTimeOfFlight_);

    // Compute (constant) mass rate.
    // double massRate = - maximumThrust_ /
    //         ( specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
    //
    // // Compute time during which the engine was switched on.
    // double engineSwitchedOnDuration = ( massAtTimeOfFlight_ - initialSpacecraftMass_ ) / massRate;
    //
    // std::shared_ptr< numerical_quadrature::QuadratureSettings< double > > quadratureSettings_ =
    //         std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( 0.0, 16 );
    //
    // // Thrust acceleration function to use quadrature.
    // // Define thrust acceleration as a function of time (to be integrated to compute the associated deltaV).
    // std::function< double( const double ) > thrustAcceleration = [ = ] ( const double currentTime ){
    //
    //     // Compute current mass.
    //     double currentMass = initialSpacecraftMass_ + massRate * currentTime;
    //
    //     // Compute and return current thrust acceleration.
    //     double currentThrustAcceleration = maximumThrust_ / currentMass;
    //
    //     return currentThrustAcceleration;
    // };
    //
    // // Create numerical quadrature from quadrature settings.
    // std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
    //         numerical_quadrature::createQuadrature( thrustAcceleration, quadratureSettings_, engineSwitchedOnDuration );
    //
    // // Compute deltaV analytically.
    // double deltaV = - specificImpulse_ * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION
    //         * std::log( 1.0 + massRate / initialSpacecraftMass_ * engineSwitchedOnDuration );
}



} // namespace low_thrust_trajectories
} // namespace tudat
