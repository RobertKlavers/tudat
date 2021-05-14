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


#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethod.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridOptimisationSetup.h"
#include "pagmo/problems/unconstrain.hpp"
#include "pagmo/algorithms/compass_search.hpp"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"
#include "Tudat/Astrodynamics/Propulsion/costateBasedThrustGuidance.h"


namespace tudat
{
namespace low_thrust_trajectories
{

HybridMethodModel HybridMethod::getModelForDecisionVector(std::vector<double> designVector) {

    // TODO Extract to common method for both here and in pagmo UDP
    Eigen::VectorXd initialCostates = Eigen::VectorXd::Zero( 6 );
    Eigen::VectorXd finalCostates = Eigen::VectorXd::Zero( 6 );

    // Check consistency of the size of the design variables vector.
    if ( designVector.size( ) != 13 )
    {
        throw std::runtime_error( "Error, size of the design variables vector unconsistent with initial and final "
                                  "MEE costates sizes." );
    }

    double timeOfFlight = designVector[0];

    for ( unsigned int i = 0 ; i < 6 ; i++ )
    {
        initialCostates( i ) = designVector[ i+1 ];
        finalCostates( i ) = designVector[ i+1 + 6 ];
    }

    // Make sure to reset the initial Mass
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass(initialMass_);

    // Create hybrid method leg.
    HybridMethodModel hybridMethodModel = HybridMethodModel(
            stateAtDeparture_, stateAtArrival_, initialCostates, finalCostates, maximumThrust_,
            specificImpulse_, timeOfFlight, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings_, hybridOptimisationSettings_ );
    return hybridMethodModel;
}

//! Perform optimisation.
std::pair< std::vector< double >, std::vector< double > > HybridMethod::performOptimisation( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 42 );

    std::string outputDirectory = "/home/robert/tud/thesis/analysis/data/hybridResults/";

    // Create object to compute the problem fitness
    problem prob{
            HybridMethodProblem(stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulse_, timeOfFlight_,
                                bodyMap_,
                                bodyToPropagate_, centralBody_, integratorSettings_, hybridOptimisationSettings_,
                                initialGuessCostates_, initialAndFinalMEEcostatesBounds_,
                                optimisationSettings_->relativeToleranceConstraints_)};

    std::vector<double> constraintsTolerance;
    for ( unsigned int i = 0 ; i < ( prob.get_nec() + prob.get_nic() ) ; i++ )
    {
        constraintsTolerance.push_back( 1.0e-3 );
    }
    prob.set_c_tol( constraintsTolerance );

    algorithm algo = optimisationSettings_->optimisationAlgorithm_;

    unsigned long long populationSize = optimisationSettings_->numberOfIndividualsPerPopulation_;

    island island{ algo, prob, populationSize };

    std::map<int, Eigen::VectorXd> hybridFitnessResults;
    std::map<int, Eigen::Vector6d> hybridErrorResults;

    // // Evolve for N generations
    for( int gen = 0 ; gen < optimisationSettings_->numberOfGenerations_ ; gen++ )
    {
        // Evolve the current population
        island.evolve( );
        while( island.status( ) != pagmo::evolve_status::idle &&
               island.status( ) != pagmo::evolve_status::idle_error )
        {
            island.wait( );
        }
        // island.wait_check( ); // Raises errors
        // Get the best individual of this generation
        std::vector<double> championDesignVector = island.get_population().champion_x();

        // Make sure to reset the initial Mass
        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass(initialMass_);

        HybridMethodModel currentHybridMethodModel = getModelForDecisionVector(championDesignVector);

        std::pair<Eigen::VectorXd, Eigen::Vector6d> fitnessResults = currentHybridMethodModel.calculateFitness();
        std::pair<std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd >> currentBestTrajectory = currentHybridMethodModel.getTrajectoryOutput();

        hybridFitnessResults[gen] = fitnessResults.first;
        hybridErrorResults[gen] = fitnessResults.second;

        double rad2deg = 180.0 /  mathematical_constants::PI;

        Eigen::Vector6d error_vec;
        error_vec <<
            fitnessResults.second(0),
            fitnessResults.second(1),
            fitnessResults.second(2) * rad2deg,
            fitnessResults.second(3) * rad2deg,
            fitnessResults.second(4) * rad2deg,
            fitnessResults.second(5) * rad2deg;

        input_output::writeDataMapToTextFile(currentBestTrajectory.first,
                                             "HybridCurrentBestTrajectory.dat",
                                             outputDirectory,
                                             "",
                                             std::numeric_limits<double>::digits10, std::numeric_limits<double>::digits10,
                                             ",");
        input_output::writeDataMapToTextFile(currentBestTrajectory.second,
                                             "HybridCurrentBestDependentVariables.dat",
                                             outputDirectory,
                                             "",
                                             std::numeric_limits<double>::digits10, std::numeric_limits<double>::digits10,
                                             ",");

        // Print results
        std::cout << "\n==== gen: " << gen << ", f: " << island.get_population().champion_f()[0] << " ====" << std::endl;
        std::cout << "  best: [" << (championDesignVector[0]/physical_constants::JULIAN_DAY) << "], ";
        for (int j = 1; j < 7; j++) {
            std::cout << "[" << championDesignVector[j] << ", " << championDesignVector[j + 6] << "], ";
        }
        std::cout << "(" << championDesignVector.size() << ")" << std::endl;
        std::cout << "  eps: [" << fitnessResults.first.transpose() << "]" << std::endl;
        std::cout << "  err: [" << error_vec.transpose() << "]" <<std::endl;
    }

    championFitness_ = island.get_population().champion_f();
    championDesignVariables_ = island.get_population().champion_x();

    HybridMethodModel championHybridMethodModel = getModelForDecisionVector(championDesignVariables_);
    std::pair<std::map< double, Eigen::VectorXd >, std::map< double, Eigen::VectorXd >> hybridChampionTrajectory = championHybridMethodModel.getTrajectoryOutput();

    // Temporary stuff to directly store the dependent variable history (hopefully containing thrust acceleration profile)
    std::cout << "Exporting Optimization Results" << std::endl;



    input_output::writeDataMapToTextFile(hybridFitnessResults,
                                         "HybridFitnessResults.dat",
                                         outputDirectory,
                                         "",
                                         std::numeric_limits<double>::digits10, std::numeric_limits<double>::digits10,
                                         ",");
    input_output::writeDataMapToTextFile(hybridErrorResults,
                                         "HybridErrorResults.dat",
                                         outputDirectory,
                                         "",
                                         std::numeric_limits<double>::digits10, std::numeric_limits<double>::digits10,
                                         ",");
    input_output::writeDataMapToTextFile(hybridChampionTrajectory.first,
                                         "HybridChampionTrajectory.dat",
                                         outputDirectory,
                                         "",
                                         std::numeric_limits<double>::digits10, std::numeric_limits<double>::digits10,
                                         ",");
    input_output::writeDataMapToTextFile(hybridChampionTrajectory.second,
                                         "HybridChampionDependentVariableHistory.dat",
                                         outputDirectory,
                                         "",
                                         std::numeric_limits<double>::digits10, std::numeric_limits<double>::digits10,
                                         ",");

    std::pair< std::vector< double >, std::vector< double > > output;
    output.first = championFitness_;
    output.second = championDesignVariables_;

    return output;

}


//! Compute current thrust vector.
Eigen::Vector3d HybridMethod::computeCurrentThrustForce(
        double time,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    std::cout<<"computeCurrentThrustForce"<<std::endl;
    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings = hybridMethodModel_->getMEEcostatesBasedThrustMagnitudeSettings( );

    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = simulation_setup::getBodyFixedThrustDirection(
                thrustMagnitudeSettings, bodyMap_, bodyToPropagate_ );

    std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction = [ = ] ( )
    {
        return computeCurrentStateVector( time );
    };

    std::function< Eigen::Vector6d( ) > centralBodyStateFunction = [ = ] ( )
    {
        return Eigen::Vector6d::Zero( );
    };

    std::function< double( ) > centralBodyGravitationalParameterFunction = [ = ]( )
    {
        return bodyMap_[ centralBody_ ]->getGravityFieldModel( )->getGravitationalParameter( );
    };

    std::function< double( ) > thrustingBodyMassFunction = std::bind( &simulation_setup::Body::getBodyMass, bodyMap_.at( bodyToPropagate_ ) );


    propulsion::MeeCostatesBangBangThrustMagnitudeWrapper thrustMagnitudeWrapper = propulsion::MeeCostatesBangBangThrustMagnitudeWrapper(
                thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                hybridMethodModel_->getCostatesFunction_( ), maximumThrust_, specificImpulseFunction, thrustingBodyMassFunction);

    propulsion::MeeCostateBasedThrustGuidance thrustGuidance = propulsion::MeeCostateBasedThrustGuidance(
                thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                hybridMethodModel_->getCostatesFunction_( ), bodyFixedThrustDirection );

    thrustGuidance.updateCalculator( time );
    thrustMagnitudeWrapper.update( time );

    return thrustMagnitudeWrapper.getCurrentThrustMagnitude( ) * thrustGuidance.getCurrentForceDirectionInPropagationFrame( );
}



//! Return thrust profile.
void HybridMethod::getThrustForceProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    std::cout<<"getThrustForceProfile"<<std::endl;
    thrustProfile.clear( );

    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > thrustMagnitudeSettings = hybridMethodModel_->getMEEcostatesBasedThrustMagnitudeSettings( );

    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = simulation_setup::getBodyFixedThrustDirection(
                thrustMagnitudeSettings, bodyMap_, bodyToPropagate_ );

    std::map< double, Eigen::Vector6d > trajectory;
    getTrajectory( epochsVector, trajectory );

    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a hybrid method trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector6d currentStateVector = trajectory[ epochsVector[ i ] ];

        std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction = [ = ] ( )
        {
            return currentStateVector;
        };

        std::function< Eigen::Vector6d( ) > centralBodyStateFunction = [ = ] ( )
        {
            return Eigen::Vector6d::Zero( );
        };

        std::function< double( ) > centralBodyGravitationalParameterFunction = [ = ]( )
        {
            return bodyMap_[ centralBody_ ]->getGravityFieldModel( )->getGravitationalParameter( );
        };

        std::function< double( ) > thrustingBodyMassFunction = std::bind( &simulation_setup::Body::getBodyMass, bodyMap_.at( bodyToPropagate_ ) );


        propulsion::MeeCostatesBangBangThrustMagnitudeWrapper thrustMagnitudeWrapper = propulsion::MeeCostatesBangBangThrustMagnitudeWrapper(
                    thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                    hybridMethodModel_->getCostatesFunction_( ), maximumThrust_, specificImpulseFunction, thrustingBodyMassFunction);

        propulsion::MeeCostateBasedThrustGuidance thrustGuidance = propulsion::MeeCostateBasedThrustGuidance(
                    thrustingBodyStateFunction, centralBodyStateFunction, centralBodyGravitationalParameterFunction,
                    hybridMethodModel_->getCostatesFunction_( ), bodyFixedThrustDirection );

        thrustGuidance.updateCalculator( epochsVector[ i ] );
        thrustMagnitudeWrapper.update( epochsVector[ i ] );

        thrustProfile[ epochsVector[ i ] ] = thrustMagnitudeWrapper.getCurrentThrustMagnitude( ) * thrustGuidance.getCurrentForceDirectionInPropagationFrame( );

    }
}

//! Compute magnitude thrust acceleration.
double HybridMethod::computeCurrentThrustAccelerationMagnitude(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    std::cout<<"computeCurrentThrustAccelerationMagnitude"<<std::endl;
    double currentMass = computeCurrentMass( currentTime, specificImpulseFunction, integratorSettings );
    Eigen::Vector3d currentThrustVector = computeCurrentThrustForce( currentTime, specificImpulseFunction, integratorSettings_ );

    return currentThrustVector.norm( ) / currentMass;

}


//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d HybridMethod::computeCurrentThrustAccelerationDirection(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    std::cout<<"computeCurrentThrustAccelerationDirection"<<std::endl;
    Eigen::Vector3d currentThrustVector = computeCurrentThrustForce( currentTime, specificImpulseFunction, integratorSettings );

    Eigen::Vector3d thrustAcceleration = currentThrustVector.normalized( );

    return thrustAcceleration.normalized( );
}



//! Return thrust acceleration profile.
void HybridMethod::getThrustAccelerationProfile(
        std::vector< double >& epochsVector,
        std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
        std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    std::cout<<"getThrustAccelerationProfile"<<std::endl;
    thrustAccelerationProfile.clear();

    std::map< double, Eigen::VectorXd > thrustProfile;
    getThrustForceProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );

    std::map< double, Eigen::VectorXd > massProfile;
    getMassProfile( epochsVector, massProfile, specificImpulseFunction, integratorSettings );

    for ( unsigned int i = 0 ; i < epochsVector.size() ; i++ )
    {
        if ( ( i > 0 ) && ( epochsVector[ i ] < epochsVector[ i - 1 ] ) )
        {
            throw std::runtime_error( "Error when retrieving the thrust profile of a hybrid method trajectory, "
                                      "epochs are not provided in increasing order." );
        }

        Eigen::Vector3d currentThrustVector = thrustProfile[ epochsVector[ i ] ];
        double currentMass = massProfile[ epochsVector[ i ] ][ 0 ];

        Eigen::Vector3d currentThrustAccelerationVector = currentThrustVector / currentMass;

        thrustAccelerationProfile[ epochsVector[ i ] ] = currentThrustAccelerationVector;

    }
}


//! Compute current cartesian state.
Eigen::Vector6d HybridMethod::computeCurrentStateVector( const double currentTime )
{
    std::cout<<"computeCurrentStateVector"<<std::endl;
    Eigen::Vector6d stateVector;
    if ( currentTime == 0.0 )
    {
        stateVector = stateAtDeparture_;
    }
    else
    {
        stateVector = hybridMethodModel_->propagateTrajectory( 0.0, currentTime, stateAtDeparture_, initialMass_).first;
    }

    return stateVector;

}


//! Retrieve acceleration map (thrust and central gravity accelerations).
basic_astrodynamics::AccelerationMap HybridMethod::retrieveLowThrustAccelerationMap(
        const simulation_setup::NamedBodyMap& bodyMapTest,
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const std::function< double ( const double ) > specificImpulseFunction,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    std::cout<<"retrieveLowThrustAccelerationMap"<<std::endl;
    basic_astrodynamics::AccelerationMap hybridMethodAccelerationMap = hybridMethodModel_->getLowThrustTrajectoryAccelerationMap( );
    return hybridMethodAccelerationMap;
}


//! Define appropriate translational state propagator settings for the full propagation.
std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >
HybridMethod::createLowThrustTranslationalStatePropagatorSettings(
        const std::string& bodyToPropagate,
        const std::string& centralBody,
        const basic_astrodynamics::AccelerationMap& accelerationModelMap,
        const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave )
{
    std::cout<<"createLowThrustTranslationalStatePropagatorSettings"<<std::endl;
    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
            std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight_, true );

    // Compute state vector at half of the time of flight.
    Eigen::Vector6d stateAtHalfOfTimeOfFlight = computeCurrentStateVector( timeOfFlight_ / 2.0 );

    // Define translational state propagator settings.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
              std::vector< std::string >{ bodyToPropagate_ }, stateAtHalfOfTimeOfFlight,
              terminationConditions.first, propagators::gauss_modified_equinoctial, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
            ( std::vector< std::string >{ centralBody_ }, accelerationModelMap,
              std::vector< std::string >{ bodyToPropagate_ }, stateAtHalfOfTimeOfFlight,
              terminationConditions.second, propagators::gauss_modified_equinoctial, dependentVariablesToSave );

    return translationalStatePropagatorSettings;
}


} // namespace low_thrust_trajectories
} // namespace tudat
