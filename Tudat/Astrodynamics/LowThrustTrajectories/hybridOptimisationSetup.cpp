/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethodModel.h"

namespace tudat
{
namespace low_thrust_trajectories
{

HybridMethodProblem::HybridMethodProblem(
        const Eigen::Vector6d &stateAtDeparture,
        const Eigen::Vector6d &stateAtArrival,
        const double maximumThrust,
        const double specificImpulse,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap bodyMap,
        const std::string bodyToPropagate,
        const std::string centralBody,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::pair< std::vector< double >, double > initialGuessThrustModel,
        const std::pair< double, double > initialAndFinalMEEcostatesBounds,
        const double relativeToleranceConstraints ) :
    stateAtDeparture_( stateAtDeparture ),
    stateAtArrival_( stateAtArrival ),
    maximumThrust_( maximumThrust ),
    specificImpulse_( specificImpulse ),
    timeOfFlight_( timeOfFlight ),
    bodyMap_( bodyMap ),
    bodyToPropagate_( bodyToPropagate ),
    centralBody_( centralBody ),
    integratorSettings_( integratorSettings ),
    initialGuessThrustModel_( initialGuessThrustModel ),
    initialAndFinalMEEcostatesBounds_( initialAndFinalMEEcostatesBounds ),
    relativeToleranceConstraints_( relativeToleranceConstraints )
{
    initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

    // Retrieve initial guess.
    guessInitialAndFinalCostates_ = initialGuessThrustModel_.first;
    relativeMarginWrtInitialGuess_ = initialGuessThrustModel_.second;

    if ( ( relativeMarginWrtInitialGuess_ != TUDAT_NAN ) && ( ( relativeMarginWrtInitialGuess_ < 0.0) || ( relativeMarginWrtInitialGuess_ > 1.0 ) ) )
    {
        std::string errorMessage = "Error when using initial guess to create hybrid method optimisation problem, the relative margin w.r.t. the "
                                   "initial guess should be constrained between [0,1], the current value is "
                + std::to_string( relativeMarginWrtInitialGuess_ ) + ".";

        throw std::runtime_error( errorMessage );
    }
}


//! Descriptive name of the problem
std::string HybridMethodProblem::get_name() const {
    return "Hybrid method to compute a low-thrust trajectory";
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > HybridMethodProblem::get_bounds() const {

    // Define lower bounds.
    std::vector< double > lowerBounds;

    // Define upper bounds.
    std::vector< double > upperBounds;

    if ( guessInitialAndFinalCostates_.size( ) != 0 )
    {
        if ( guessInitialAndFinalCostates_.size() != 12 )
        {
            throw std::runtime_error( "Error when providing an initial guess for hybrid method, size of the vector unconsistent"
                                      "with the expected 5 initial and 5 final MEE costate values." );
        }
        else
        {
            for ( int i = 0 ; i < 12 ; i++ )
            {
                double lowerBoundsFromInitialGuess = ( 1.0 - relativeMarginWrtInitialGuess_ ) * guessInitialAndFinalCostates_[ i ];
                lowerBounds.push_back( lowerBoundsFromInitialGuess );
            }

            for ( int i = 0 ; i < 12 ; i++ )
            {
                double upperBoundsFromInitialGuess = ( 1.0 + relativeMarginWrtInitialGuess_ ) * guessInitialAndFinalCostates_[ i ];
                upperBounds.push_back( upperBoundsFromInitialGuess );
            }
        }
    }
    else
    {
        for ( int i = 0 ; i < 12 ; i ++ )
        {
            lowerBounds.push_back( initialAndFinalMEEcostatesBounds_.first );
            upperBounds.push_back( initialAndFinalMEEcostatesBounds_.second );
        }
    }

    return { lowerBounds, upperBounds };
}

//! Fitness function.
std::vector< double > HybridMethodProblem::fitness( const std::vector< double > &designVariables ) const {

    // Re-initialise mass of the spacecraft.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Transform vector of design variables into 3D vector of throttles.
    Eigen::VectorXd initialCostates = Eigen::VectorXd::Zero( 6 );
    Eigen::VectorXd finalCostates = Eigen::VectorXd::Zero( 6 );

    // Check consistency of the size of the design variables vector.
    if ( designVariables.size( ) != 12 )
    {
        throw std::runtime_error( "Error, size of the design variables vector unconsistent with initial and final "
                                  "MEE costates sizes." );
    }


    for ( unsigned int i = 0 ; i < 6 ; i++ )
    {
        initialCostates( i ) = designVariables[ i ];
        finalCostates( i ) = designVariables[ i + 6 ];
    }

    std::vector< double > fitness;

    // Create hybrid method leg.
    low_thrust_trajectories::HybridMethodModel currentLeg = low_thrust_trajectories::HybridMethodModel(
                stateAtDeparture_, stateAtArrival_, initialCostates, finalCostates, maximumThrust_,
                specificImpulse_, timeOfFlight_, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings_ );

    // Propagate until time of flight is reached.
    Eigen::Vector6d finalPropagatedState = currentLeg.propagateTrajectory( );

    // Convert final propagated state to MEE.
    Eigen::Vector6d finalPropagatedKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
                finalPropagatedState, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );

    // Convert targeted final state to MEE.
    Eigen::Vector6d finalTargetedKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
                stateAtArrival_, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );

    Eigen::Vector6d error = (finalPropagatedKeplerianState - finalTargetedKeplerianState).cwiseAbs();

    // Compute auxiliary variables.
    Eigen::Vector6d c;
    Eigen::Vector6d r;
    Eigen::Vector6d epsilon;
    Eigen::Vector6d epsilon_lower;
    Eigen::Vector6d epsilon_upper;

    // TODO Read weights and constraints from input configuration
    Eigen::Vector6d constraint_weights;
    constraint_weights << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
    double weightMass = 1000.0;
    double weightTOF = 0.0;

    epsilon_lower = Eigen::Vector6d::Zero();
    epsilon_upper <<
        100.0e3,
        0.01,
        0.1 * mathematical_constants::PI / 180.0,
        1.0 * mathematical_constants::PI / 180.0,
        1.0 * mathematical_constants::PI / 180.0,
        0.0 * mathematical_constants::PI / 180.0;

    // 5 for no longitude targeting?
    int finalErrorParameterCount = 3;
    for ( int i = 0 ; i < finalErrorParameterCount ; i++ )
    {
        c[ i ] = 1.0 / (epsilon_upper(i) - epsilon_lower(i));
        r[ i ] = 1.0 - epsilon_upper(i) * c[ i ];
        epsilon[ i ] = error[ i ] * c[ i ] + r[ i ];
    }
    // Determine scaled epsilon at TOF
    double epsilon_final = 0.0;
    for ( int i = 0 ; i < finalErrorParameterCount ; i++) {
        epsilon_final += constraint_weights[i] * epsilon[i] * epsilon[i];
    }

    //TODO ADD AS DEBUG OUTPUT
    // std::cout << "err: " << error.transpose() << std::endl;
    double finalMass = currentLeg.getMassAtTimeOfFlight();

    // AOF From Jimenez: F = sum(orbit error) + W_t*t_f + W_m*(1 - m_f/m_0)
    double optimisationObjectiveHostep = epsilon_final + weightTOF * timeOfFlight_ + weightMass*(1 - finalMass/initialSpacecraftMass_);

    // std::cout << "f: "<< optimisationObjectiveHostep << ", eps: " << epsilon.transpose() << std::endl;

    // Output of the fitness function..
    fitness.push_back( optimisationObjectiveHostep );

    return fitness;
}

} // namespace low_thrust_trajectories
} // namespace tudat


