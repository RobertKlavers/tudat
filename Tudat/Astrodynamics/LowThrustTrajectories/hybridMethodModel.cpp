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
#include <Tudat/SimulationSetup/hybridOptimisationSettings.h>
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethodModel.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
int epochCount = 0;
namespace tudat
{
namespace low_thrust_trajectories
{

using namespace orbital_element_conversions;

std::shared_ptr< simulation_setup::AccelerationSettings > HybridMethodModel::getTangentialThrustAccelerationSettings() {
    std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< simulation_setup::ThrustDirectionFromStateGuidanceSettings >( centralBody_, true, false );
    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< simulation_setup::ConstantThrustMagnitudeSettings >( maximumThrust_, specificImpulse_ );

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsOfVehicle;
    return std::make_shared< simulation_setup::ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings );
}

std::shared_ptr< simulation_setup::AccelerationSettings > HybridMethodModel::getRadialThrustAccelerationSettings() {
    std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< simulation_setup::ThrustDirectionFromStateGuidanceSettings >( centralBody_, false, false );
    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< simulation_setup::ConstantThrustMagnitudeSettings >( maximumThrust_, specificImpulse_ );

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsOfVehicle;
    return std::make_shared< simulation_setup::ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings );
}

std::shared_ptr< simulation_setup::AccelerationSettings > HybridMethodModel::getOutOfPlaneThrustAccelerationSettings() {
    // getForceDirectionOutOfPlane

    std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< simulation_setup::ThrustDirectionFromStateGuidanceSettings >( centralBody_, false, true, true );

    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings =
            std::make_shared< simulation_setup::ConstantThrustMagnitudeSettings >( maximumThrust_, specificImpulse_ );

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > accelerationsOfVehicle;
    return std::make_shared< simulation_setup::ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings );
}


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

    switch(hybridOptimisationSettings_->propagationType_) {
        case simulation_setup::tangential: {
            bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getTangentialThrustAccelerationSettings( ) );
            break;
        }
        case simulation_setup::radial: {
            bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getRadialThrustAccelerationSettings( ) );
            break;
        }
        case simulation_setup::outofplane: {
            bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getOutOfPlaneThrustAccelerationSettings( ) );
            break;
        }
        case simulation_setup::costates: {
            bodyToPropagateAccelerations[ bodyToPropagate_ ].push_back( getMEEcostatesBasedThrustAccelerationSettings( ) );
        }
    }

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ bodyToPropagate_ ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap_, accelerationMap, std::vector< std::string >{ bodyToPropagate_ }, std::vector< std::string >{ centralBody_ } );

    return accelerationModelMap;
}

bool customTerminationFunction(double time, int numberOfEpochs) {
    epochCount++;
    return (epochCount >= numberOfEpochs);
}

propagators::SingleArcDynamicsSimulator<> HybridMethodModel::getDynamicsSimulator(
        double initialTime,
        double finalTime,
        Eigen::Vector6d initialState,
        double initialMass,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        bool withDependent,
        bool useOA,
        int numberOfSteps
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

    // Termination Settings set up
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings;

    // Terminate when periapsis is too low
    std::shared_ptr< propagators::SingleDependentVariableSaveSettings > altitudeTerminationVariable =
            std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::periapsis_altitude_dependent_variable, bodyToPropagate_, centralBody_ );

    // Terminate when eccentricity is too large
    std::shared_ptr< propagators::SingleDependentVariableSaveSettings > eccentricityTerminationVariable =
            std::make_shared<propagators::SingleDependentVariableSaveSettings>(
                    propagators::keplerian_state_dependent_variable, bodyToPropagate_, centralBody_, 1);

    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettingsList;

    // TODO: Configure termination boundaries 'magic' values
    if (!useOA) {
        terminationSettingsList.push_back(std::make_shared< propagators::PropagationTimeTerminationSettings >( finalTime, false ));
    }
    terminationSettingsList.push_back(std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(altitudeTerminationVariable, 150.0e3, 1, false));
    terminationSettingsList.push_back(std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(eccentricityTerminationVariable, 0.9, 0, false));

    if (useOA) {
        //TODO: From Configuration
        epochCount = 0;
        std::shared_ptr< propagators::PropagationTerminationSettings > customTerminationSettings =
                std::make_shared< propagators::PropagationCustomTerminationSettings >(
                        std::bind( &customTerminationFunction, std::placeholders::_1, numberOfSteps ) );
        terminationSettingsList.push_back(customTerminationSettings);
    }

    terminationSettings = std::make_shared< propagators::PropagationHybridTerminationSettings >(
            terminationSettingsList, true );

    // Dependent variables set up
    std::vector<std::shared_ptr<propagators::SingleDependentVariableSaveSettings> > dependentVariablesList;

    dependentVariablesList.push_back(std::make_shared<propagators::SingleDependentVariableSaveSettings>(
            propagators::keplerian_state_dependent_variable, bodyToPropagate_, centralBody_));
    dependentVariablesList.push_back(std::make_shared<propagators::SingleAccelerationDependentVariableSaveSettings>(
            basic_astrodynamics::thrust_acceleration, bodyToPropagate_, bodyToPropagate_, 0));
    dependentVariablesList.push_back(std::make_shared<propagators::SingleDependentVariableSaveSettings>(
            propagators::lvlh_to_inertial_frame_rotation_dependent_variable, bodyToPropagate_, centralBody_));

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
    Eigen::VectorXd propagationResult = propagateTrajectory( 0.0, timeOfFlight_, stateAtDeparture_, initialSpacecraftMass_ ).first;
    massAtTimeOfFlight_ = propagationResult[ 6 ];
    return propagationResult.segment( 0, 6 );
}

std::map<double, Eigen::VectorXd> HybridMethodModel::propagateTrajectoryBenchmark(double stepSize) {
    // Set up specific settings for this OA Arc
    std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings =
            std::make_shared<numerical_integrators::IntegratorSettings<double> >
                    (numerical_integrators::rungeKutta4, 0.0, stepSize);

    // Re-initialise integrator settings.
    bodyMap_[bodyToPropagate_]->setConstantBodyMass(initialSpacecraftMass_);
    basic_astrodynamics::AccelerationMap accelerationModelMap = getLowThrustTrajectoryAccelerationMap();

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModel;
    massRateModel[ bodyToPropagate_ ] = createMassRateModel( bodyToPropagate_, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                             bodyMap_, accelerationModelMap );

    // Termination Settings set up
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings;
    terminationSettings = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, false );

    // Define propagator settings.
    std::shared_ptr<propagators::TranslationalStatePropagatorSettings<double> > translationalStatePropagatorSettings =
            std::make_shared<propagators::TranslationalStatePropagatorSettings<double> >(
                    std::vector<std::string>{centralBody_}, accelerationModelMap,
                    std::vector<std::string>{bodyToPropagate_},
                    stateAtDeparture_, terminationSettings, propagators::cowell);

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
                                                                                terminationSettings);

    // Propagate the trajectory.
    propagators::SingleArcDynamicsSimulator<> dynamicsSimulator(bodyMap_, integratorSettings, propagatorSettings);

    return dynamicsSimulator.getEquationsOfMotionNumericalSolution();
}

std::tuple<double, Eigen::VectorXd, Eigen::VectorXd> HybridMethodModel::getStateIncrease(double initialTime, Eigen::Vector6d initialState, double initialMass, int numberOfSteps, double numberOfRevolutionsToPropagate) {
    // std::cout << "t0: " << initialTime << "\nx0: " << initialState.transpose() << "\nm0: " << initialMass << std::endl;
    // Set up initial conditions for the first OA Arc
    double stepSize = 2.0 * mathematical_constants::PI / numberOfSteps;
    // Set up specific settings for this OA Arc
    std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings =
            std::make_shared<numerical_integrators::IntegratorSettings<double> >
                    (numerical_integrators::rungeKutta4, initialTime, stepSize);

    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator = getDynamicsSimulator(initialTime, timeOfFlight_, initialState, initialMass, integratorSettings, false, true, numberOfSteps);

    // Check if the propagation was terminated for the expected single orbital evolution termination condition, stop
    // execution completely
    propagators::PropagationTerminationReason terminationReason = dynamicsSimulator.getPropagationTerminationReason()->getPropagationTerminationReason();
    std::shared_ptr< propagators::HybridPropagationTerminationCondition > hybridPropagationTerminationCondition =
            std::dynamic_pointer_cast< propagators::HybridPropagationTerminationCondition >( dynamicsSimulator.getPropagationTerminationCondition( ) );

    // Expected order is: 0) time termination, 1) altitude termination, 2) eccentricity termination, 3) custom termination
    std::vector<bool> terminationConditionsMet = hybridPropagationTerminationCondition->getIsConditionMetWhenStopping();

    bool terminatedOnSingleRevolution = (
            terminationReason == propagators::termination_condition_reached &&
            terminationConditionsMet[2]);

    if (!terminatedOnSingleRevolution) {
        throw std::runtime_error( "Error, propagation not terminated by OA" );
    }

    // Integrate single OA Arc, not that MEE are expected
    std::map< double, Eigen::VectorXd > rawNumericalSolution = dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );


    // Keep track of the state derivatives and time epochs
    std::vector<Eigen::VectorXd> stateDerivatives;
    std::vector<double> timeEpochs;

    for( std::pair<double, Eigen::VectorXd> element : rawNumericalSolution ) {
        Eigen::VectorXd stateDerivative = dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(element.first, element.second);
        stateDerivatives.emplace_back(stateDerivative);
        timeEpochs.emplace_back(element.first);
    }

    // Use Trapezoidal Quadrature Integrator (in time) to integrate states
    tudat::numerical_quadrature::TrapezoidNumericalQuadrature< double, Eigen::VectorXd > integrator(timeEpochs, stateDerivatives);
    Eigen::VectorXd computedIntegralTrapezoid = integrator.getQuadrature();

    // Calculate the orbital period from the calculated time steps
    double orbitalPeriod = (timeEpochs.back() - timeEpochs.front());

    // Averages for this revolution
    // Eigen::VectorXd averageStateProgression = computedIntegralTrapezoid / orbitalPeriod;
    // Eigen::VectorXd averageStateProgression = computedIntegralTrapezoid;

    // double finalOAtime = initialTime + numberOfRevolutionsToPropagate * orbitalPeriod;

    return {orbitalPeriod, rawNumericalSolution.cbegin()->second, computedIntegralTrapezoid/orbitalPeriod};
}


std::map<double, Eigen::Vector6d> HybridMethodModel::propagateTrajectoryOA(double numberOfRevs, int numberOfSteps ) {
    double initialArcTime = 0.0;
    double initialArcMass = initialSpacecraftMass_;
    Eigen::Vector6d initialArcStateCartesian = stateAtDeparture_;

    // Retrieve gravitational parameter for convenience
    double gravitationalParameter = bodyMap_[centralBody_]->getGravityFieldModel()->getGravitationalParameter();

    // Store trajectory results
    std::map< double, Eigen::Vector6d > resultOATrajectory;
    resultOATrajectory[initialArcTime] = stateAtDeparture_;

    // Minimum number of OA Arcs needed to reach TOF
    // int numberOfOAarcs = floor(timeOfFlight_ / rawIntermediateState);

    bool finalTimeReached = false;

    while (!finalTimeReached) {

        double intermediateOrbitalPeriod;
        Eigen::VectorXd rawInitialArcState;
        Eigen::VectorXd rawIntermediateStateIncrease;
        try {
            std::tie(intermediateOrbitalPeriod, rawInitialArcState, rawIntermediateStateIncrease) = getStateIncrease(initialArcTime, initialArcStateCartesian, initialArcMass, numberOfSteps, numberOfRevs);
        } catch (std::runtime_error& error) {
            break;
        }

        double numberOfRevolutionsToPropagate = numberOfRevs;
        bool propagationReachesTimeOfFlight = initialArcTime + numberOfRevs * intermediateOrbitalPeriod > timeOfFlight_;
        double timeStepOA = numberOfRevolutionsToPropagate * intermediateOrbitalPeriod;

        if (propagationReachesTimeOfFlight) {
            timeStepOA = timeOfFlight_ - initialArcTime;
            finalTimeReached = true;
        }

        double intermediateTime = initialArcTime + timeStepOA;
        Eigen::VectorXd rawIntermediateState = rawInitialArcState + timeStepOA*rawIntermediateStateIncrease;
        Eigen::Vector6d intermediateOAMEEState = rawIntermediateState.segment(0, 6);
        Eigen::Vector6d intermediateOACartesianState = orbital_element_conversions::convertModifiedEquinoctialToCartesianElements(intermediateOAMEEState, gravitationalParameter, false);
        double intermediateOASpacecraftMass = rawIntermediateState[6];

        double finalOrbitalPeriod;
        Eigen::VectorXd rawFinalArcState;
        Eigen::VectorXd rawFinalStateIncrease;
        try {
            std::tie(finalOrbitalPeriod, rawFinalArcState, rawFinalStateIncrease) = getStateIncrease(intermediateTime, intermediateOACartesianState, intermediateOASpacecraftMass, numberOfSteps, numberOfRevs);
        } catch (std::runtime_error& error) {
            break;
        }

        Eigen::VectorXd rawFinalOAArcState = rawInitialArcState + (timeStepOA / 2) * (rawIntermediateStateIncrease + rawFinalStateIncrease);

        Eigen::Vector6d finalOAMEEState = rawFinalOAArcState.segment(0,6);
        Eigen::Vector6d finalOAArcState = orbital_element_conversions::convertModifiedEquinoctialToCartesianElements(finalOAMEEState, gravitationalParameter, false);
        double finalOAArcMass = rawFinalOAArcState[6];

        initialArcTime = intermediateTime;
        initialArcStateCartesian = finalOAArcState;
        initialArcMass = finalOAArcMass;
        resultOATrajectory[intermediateTime] = initialArcStateCartesian;
    }
    massAtTimeOfFlight_ = initialArcMass;
    return resultOATrajectory;
}


//! Propagate the spacecraft trajectory to a given time.
std::pair<Eigen::VectorXd, Eigen::Vector6d> HybridMethodModel::propagateTrajectory( double initialTime, double finalTime, Eigen::Vector6d initialState, double initialMass)
{
    integratorSettings_->initialTime_ = initialTime;
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator = getDynamicsSimulator(initialTime, finalTime, initialState, initialMass, integratorSettings_);

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
    Eigen::VectorXd propagationOriginalResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second;

    Eigen::Vector6d computedMMEStateDerivatives;

    if ( finalTime == timeOfFlight_ )
    {
        massAtTimeOfFlight_ = propagationOriginalResult[ 6 ];
    }

    return {propagationOriginalResult, computedMMEStateDerivatives};
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
                Eigen::VectorXd propagationResult = propagateTrajectory( 0.0, currentTime, propagatedState, currentMass).first;
                propagatedState = propagationResult.segment(0, 6);
                currentMass = propagationResult[6];}
            propagatedTrajectory[ currentTime ] = propagatedState;
        }
        else
        {
            Eigen::VectorXd propagationResult = propagateTrajectory( 0.0, currentTime, propagatedState, currentMass).first;
            propagatedState = propagationResult.segment(0, 6);
            currentMass = propagationResult[6];
            propagatedTrajectory[ currentTime ] = propagatedState;
        }

    }

    bodyMap_[ centralBody_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    return propagatedTrajectory;
}

//! Utility to retrieve the integration result and all saved dependent variables
std::pair<std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd >> HybridMethodModel::getTrajectoryOutput() {
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator = getDynamicsSimulator(0.0, timeOfFlight_, stateAtDeparture_, initialSpacecraftMass_, integratorSettings_, false);
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    Eigen::VectorXd finalPropagatedState = integrationResult.rbegin( )->second;
    massAtTimeOfFlight_ = finalPropagatedState[6];

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
    if (hybridOptimisationSettings_->debug_) {
        std::cout << "  --mod.fitcalc--" << std::endl;
    }

    // Propagate until time of flight is reached.
    // Eigen::Vector6d finalPropagatedState = propagateTrajectory( );

    std::map<double, Eigen::Vector6d> finalOAResult = propagateTrajectoryOA(13, 30);
    Eigen::Vector6d finalPropagatedState = finalOAResult.rbegin()->second;

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

    // Create vector with individual contributions to objective function
    // Using AOF From (D. Jimenez, 2020): F = sum(W_j*epsilon_j^2) + W_t*t_f + W_m*(1 - m_f/m_0)
    for ( int i = 0 ; i < 6 ; i++) {
        aof_vector(i) = (hybridOptimisationSettings_->constraintWeights_[i] * epsilon[i] * epsilon[i]);
    }

    // Time of flight contribution (t_f should be in days)
    double epsTimeOfFlight = (hybridOptimisationSettings_->weightTimeOfFlight_ * timeOfFlight_/physical_constants::JULIAN_DAY);
    aof_vector(6) = epsTimeOfFlight;

    // Final mass contribution
    double finalMass = getMassAtTimeOfFlight();
    double epsFinalMass = (hybridOptimisationSettings_->weightMass_*(1 - finalMass/initialSpacecraftMass_));
    aof_vector(7) = epsFinalMass;

    // Some debugging statements we don't actually directly output the fitness
    if (hybridOptimisationSettings_->debug_) {
        double totalEps = aof_vector.sum();
        std::cout << "    cfun : [" << costatesFunction_(0.0).transpose() << "] | [" << costatesFunction_(timeOfFlight_).transpose() << "]" << std::endl;
        std::cout << "    x_f  : " << finalPropagatedState.transpose() << std::endl;
        std::cout << "    err  : " << error.transpose() << std::endl;
        std::cout << "    eps  : " << aof_vector.transpose() <<std::endl;
        std::cout << "    sum_e: " << totalEps << std::endl;
        std::cout << "  --end.model--" << std::endl;
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
