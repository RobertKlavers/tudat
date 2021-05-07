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

#ifndef TUDAT_HYBRID_METHOD_H
#define TUDAT_HYBRID_METHOD_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethodModel.h"
#include "Tudat/SimulationSetup/optimisationSettings.h"

namespace tudat
{
namespace low_thrust_trajectories
{

class HybridMethod : public low_thrust_trajectories::LowThrustLeg
{
public:

    //! Constructor.
    HybridMethod(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double centralBodyGravitationalParameter,
            const double initialSpacecraftMass,
            const double maximumThrust,
            const double specificImpulse,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap& bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings,
            const std::pair< double, double > initialAndFinalMEEcostatesBounds = std::make_pair( - 1.0e4, 1.0e4 ) ):
        LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, initialSpacecraftMass, true ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        bodyMap_( bodyMap ),
        bodyToPropagate_( bodyToPropagate ),
        centralBody_(centralBody ),
        maximumThrust_( maximumThrust ),
        specificImpulse_( specificImpulse ),
        integratorSettings_( integratorSettings ),
        optimisationSettings_( optimisationSettings ),
        initialAndFinalMEEcostatesBounds_( initialAndFinalMEEcostatesBounds )
    {

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the hybrid method.
        if ( optimisationSettings_->initialGuessThrustModel_.first.size( ) != 0 )
        {
            if ( optimisationSettings_->initialGuessThrustModel_.first.size( ) != 12 )
            {
                throw std::runtime_error( "Error when providing an initial guess for hybrid method, vector size"
                                          "unconsistent with the 5 initial and 5 final MEE costate values which are required." );
            }
            else
            {
                Eigen::VectorXd initialGuessMEEinitialAndFinalCostates( 12 ); // No entirely sure why this is needed
                initialGuessThrustModel_.first = optimisationSettings_->initialGuessThrustModel_.first;
            }
        }
        else
        {
            initialGuessThrustModel_.first = std::vector< double>( );
        }
        initialGuessThrustModel_.second = optimisationSettings_->initialGuessThrustModel_.second;

       // Perform optimisation
       std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
       championFitness_ = bestIndividual.first;
       championDesignVariables_ = bestIndividual.second;


       // Transform vector of design variables into 3D vector of throttles.
       Eigen::VectorXd initialCostates; initialCostates.resize( 6 );
       Eigen::VectorXd finalCostates; finalCostates.resize( 6 );
       for ( unsigned int i = 0 ; i < 6; i++ )
       {
           initialCostates[ i ] = championDesignVariables_[ i ];
           finalCostates[ i ] = championDesignVariables_[ i + 6 ];
       }

        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialMass_ );

        // Create Hybrid leg from the best optimisation individual.
        hybridMethodModel_ = std::make_shared< HybridMethodModel >(
                    stateAtDeparture_, stateAtArrival_, initialCostates, finalCostates, maximumThrust_, specificImpulse_, timeOfFlight_,
                    bodyMap_, bodyToPropagate_, centralBody_, integratorSettings );

    }

    //! Default destructor.
    ~HybridMethod( ) { }

    //! Convert time to independent variable.
    double convertTimeToIndependentVariable( const double time )
    {
        return time;
    }

    //! Convert independent variable to time.
    double convertIndependentVariableToTime( const double independentVariable )
    {
        return independentVariable;
    }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return hybridMethodModel_->getTotalDeltaV( );
    }

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double currentTime );

    //! Compute state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory )
    {
        propagatedTrajectory = hybridMethodModel_->propagateTrajectory( epochsVector, propagatedTrajectory );
    }


    Eigen::Vector3d computeCurrentThrustForce( double time,
                                          std::function< double ( const double ) > specificImpulseFunction,
                                          std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust profile.
    void getThrustForceProfile( std::vector< double >& epochsVector,
                           std::map< double, Eigen::VectorXd >& thrustProfile,
                           std::function< double ( const double ) > specificImpulseFunction,
                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute magnitude thrust acceleration.
    double computeCurrentThrustAccelerationMagnitude(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust acceleration profile.
    void getThrustAccelerationProfile(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Return best individual.
    std::vector< double > getBestIndividual( )
    {
        return championDesignVariables_;
    }

    //! Return fitness of best individual.
    std::vector< double > getBestIndividualFitness( )
    {
        return championFitness_;
    }

    //! Return best hybrid method leg after optimisation.
    std::shared_ptr< HybridMethodModel > getOptimalHybridMethodModel( )
    {
        return hybridMethodModel_;
    }


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
            const simulation_setup::NamedBodyMap& bodyMapTest,
            const std::string& bodyToPropagate,
            const std::string& centralBody,
            const std::function< double ( const double ) > specificImpulseFunction,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            const std::string& bodyToPropagate,
            const std::string& centralBody,
            const basic_astrodynamics::AccelerationMap& accelerationModelMap,
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


private:

    double centralBodyGravitationalParameter_;

    //! Body map object.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    double specificImpulse_;

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Optimisation settings.
    std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< std::vector< double >, double > initialGuessThrustModel_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Hybrid method leg corresponding to the best optimisation output.
    std::shared_ptr< HybridMethodModel > hybridMethodModel_;

    //! Lower and upper bounds for the initial and final MEE costates.
    std::pair< double, double > initialAndFinalMEEcostatesBounds_;


};


} // namespace low_thrust_trajectories
} // namespace tudat

#endif // TUDAT_HYBRID_METHOD_H
