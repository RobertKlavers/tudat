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
#include "Tudat/SimulationSetup/hybridOptimisationSettings.h"

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
        std::shared_ptr< simulation_setup::HybridOptimisationSettings > hybridOptimisationSettings,
        const std::pair< std::vector< double >, double > initialGuessCostates_,
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
    hybridOptimisationSettings_(hybridOptimisationSettings),
    initialGuessCostates_( initialGuessCostates_ ),
    initialAndFinalMEEcostatesBounds_( initialAndFinalMEEcostatesBounds ),
    relativeToleranceConstraints_( relativeToleranceConstraints )
{
    initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

    // Retrieve initial guess.
    guessInitialAndFinalCostates_ = initialGuessCostates_.first;
    relativeMarginWrtInitialGuess_ = initialGuessCostates_.second;

    if ( ( relativeMarginWrtInitialGuess_ != TUDAT_NAN ) && relativeMarginWrtInitialGuess_ < 0 )
    {
        std::string errorMessage = "Error when using initial guess to create hybrid method optimisation problem, the relative margin w.r.t. the "
                                   "initial guess should be larger than 0, the current value is "
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

    //TODO Add ability to set lower bound for time of flight
    lowerBounds.push_back(0.0);
    upperBounds.push_back(timeOfFlight_);

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
                double lowerBoundsFromInitialGuess = guessInitialAndFinalCostates_[ i ] - relativeMarginWrtInitialGuess_;
                lowerBounds.push_back( lowerBoundsFromInitialGuess );
            }

            for ( int i = 0 ; i < 12 ; i++ )
            {
                double upperBoundsFromInitialGuess = guessInitialAndFinalCostates_[ i ] + relativeMarginWrtInitialGuess_;;
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
    std::vector< double > fitness;

    if (hybridOptimisationSettings_->debug_) {
        std::cout << "\n--prob.fitfunc--" << std::endl;
        std::cout << "  y_vec: [";
        for (auto& n: designVariables) {
            std::cout << n << ", ";
        }
        std::cout << "]" << std::endl;
    }

    double tofDecisionVector = designVariables[0];

    // Transform vector of design variables into 3D vector of throttles.
    Eigen::VectorXd initialCostates = Eigen::VectorXd::Zero( 6 );
    Eigen::VectorXd finalCostates = Eigen::VectorXd::Zero( 6 );

    // Check consistency of the size of the design variables vector.
    if ( designVariables.size( ) != 13 )
    {
        throw std::runtime_error( "Error, size of the design variables vector unconsistent with initial and final "
                                  "MEE costates sizes." );
    }

    for ( unsigned int i = 0 ; i < 6 ; i++ ) {
        initialCostates(i) = designVariables[i + 1];
        finalCostates(i) = designVariables[i + 1 + 6];
    }

    // Re-initialise mass of the spacecraft.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Create hybrid method leg.
    HybridMethodModel currentLeg = HybridMethodModel(
            stateAtDeparture_, stateAtArrival_, initialCostates, finalCostates, maximumThrust_,
            specificImpulse_, tofDecisionVector, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings_, hybridOptimisationSettings_ );

    std::pair<Eigen::VectorXd, Eigen::Vector6d> fitnessResults = currentLeg.calculateFitness();

    double totalEpsilon = fitnessResults.first.sum();

    if (hybridOptimisationSettings_->debug_) {
        std::cout << "  f: "<< totalEpsilon << std::endl;
        std::cout << "  eps: [" << fitnessResults.first.transpose() << "]" << std::endl;
        std::cout << "  err: [" << fitnessResults.second.transpose() << "]" << std::endl;
        std::cout << "  cst: [" << initialCostates.transpose() << "] | [" << finalCostates.transpose() << "]" << std::endl;
        std::cout << "  dvc: [" << (designVariables[0]/physical_constants::JULIAN_DAY) << "], ";
        for (int j = 1; j < 7; j++) {
            std::cout << "[" << designVariables[j] << ", " << designVariables[j + 6] << "], ";
        }
        std::cout << "(" << designVariables.size() << ")" << std::endl;
        std::cout << "--end.prob--" << std::endl;
    }
    // Output of the fitness function..
    fitness.push_back( totalEpsilon );

    return fitness;
}

} // namespace low_thrust_trajectories
} // namespace tudat


