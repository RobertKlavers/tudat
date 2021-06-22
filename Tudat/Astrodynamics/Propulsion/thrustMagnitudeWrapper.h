/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_THRUSTMAGNITUDEWRAPPER_H
#define TUDAT_THRUSTMAGNITUDEWRAPPER_H

#include <memory>
#include <functional>
#include <boost/lambda/lambda.hpp>
#include <iostream>

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"

namespace tudat
{

namespace propulsion
{

//! Base class for determining the current thrust magnitude from method defined in derived class.
class ThrustMagnitudeWrapper
{
public:

    //! Constructor
    ThrustMagnitudeWrapper( ):
        currentTime_( TUDAT_NAN ){ }

    //! Destructor.
    virtual ~ThrustMagnitudeWrapper( ){ }

    //! Pure virtual function to update the thrust magnitude to the current time.
    /*!
     * Pure virtual function to update the thrust magnitude to the current time. Method of computation is to be
     * implemented in the derived class.
     * \param time Time at which the magitude is to be computed.
     */
    virtual void update( const double time ) = 0;

    //! Pure virtual function to get the current thrust magnitude.
    /*!
     * Pure virtual function to get the current thrust magnitude.
     * \return Current thrust magnitude.
     */
    virtual double getCurrentThrustMagnitude( ) = 0;

    //! Pure virtual function to get the current mass rate.
    /*!
     * Pure virtual function to get the current mass rate.
     * \return Current mass rate.
     */
    virtual double getCurrentMassRate( ) = 0;

    //! Function to reset the current time of the thrust model
    /*!
     *  Function to reset the current time of the thrust model. Function is typically used to reset the time to NaN,
     *  signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    void resetCurrentTime( const double currentTime = TUDAT_NAN )
    {
        currentTime_ = currentTime;
        resetDerivedClassCurrentTime( currentTime );
    }

    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities. This function can be redefined in
     *  derived class
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    {

    }


protected:

    //! Current time for model.
    double currentTime_;
};

//! Class for custom computations of thrust magnitude and mass rate.
/*!
 *  Class for custom computations of thrust magnitude and mass rate. Both the magnitude and specific impulse
 *  are parameterized as a function of time and provided as input, providing a flexible interface for user-defined
 *  thrust magnitude profiles.
 */
class CustomThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustMagnitudeFunction Function returning thrust as a function of time.
     * \param specificImpulseFunction Function returning specific impulse as a function of time.
     * \param isEngineOnFunction Function returning whether the function is on (returns true if so) at a given time.
     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
     * being started (empty by default)
     */
    CustomThrustMagnitudeWrapper(
            const std::function< double( const double ) > thrustMagnitudeFunction,
            const std::function< double( const double ) > specificImpulseFunction,
            const std::function< bool( const double ) > isEngineOnFunction = [ ]( const double ){ return true; } ,
            const std::function< void( const double ) > customThrustResetFunction = std::function< void( const double ) >( ) ):
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        isEngineOnFunction_( isEngineOnFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        customThrustResetFunction_( customThrustResetFunction ){ }

    //! Destructor.
    ~CustomThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            // If engine is one, update engine.
            if( isEngineOnFunction_( time ) )
            {
                currentThrustMagnitude_ = thrustMagnitudeFunction_( time );
                currentSpecificImpulse_ = specificImpulseFunction_( time );
            }
            // If engine is off, set thrust to zero.
            else
            {
                currentThrustMagnitude_ = 0.0;
                currentSpecificImpulse_ = TUDAT_NAN;
            }
            currentTime_ = time;
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustMagnitude( )
    {
        return currentThrustMagnitude_;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( )
    {
        if( currentThrustMagnitude_ != 0.0 )
        {
            return propulsion::computePropellantMassRateFromSpecificImpulse(
                        currentThrustMagnitude_, currentSpecificImpulse_ );
        }
        else
        {
            return 0.0;
        }
    }


    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    {
        if( !( customThrustResetFunction_ == nullptr ) )
        {
            customThrustResetFunction_( currentTime );\
        }
    }

private:

    //! Function returning thrust as a function of time..
    std::function< double( const double ) > thrustMagnitudeFunction_;

    //! Function returning specific impulse as a function of time.
    std::function< double( const double ) > specificImpulseFunction_;

    //! Function returning whether the function is on (returns true if so) at a given time.
    std::function< bool( const double ) > isEngineOnFunction_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

    std::function< void( const double ) > customThrustResetFunction_;

};



//! Class for bang-bang thrust magnitude from MEE co-states (optimal control theory).
/*!
 *  Class for bang-bang thrust magnitude from MEE co-states (optimal control theory).
 */
class MeeCostatesBangBangThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustMagnitudeFunction Function returning thrust as a function of time.
     * \param specificImpulseFunction Function returning specific impulse as a function of time.
     * \param isEngineOnFunction Function returning whether the function is on (returns true if so) at a given time.
     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
     * being started (empty by default)
     */
    MeeCostatesBangBangThrustMagnitudeWrapper(
            const std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction,
            const std::function< Eigen::Vector6d( ) > centralBodyStateFunction,
            const std::function< double( ) > centralBodyGravitationalParameterFunction,
            const std::function< Eigen::VectorXd( const double ) > costateFunction,
            const double thrustMagnitude,
            const std::function< double( const double ) > specificImpulseFunction,
            const std::function< double( ) > thrustingBodyMassFunction,
            const std::function< void( const double ) > customThrustResetFunction = std::function< void( const double ) >( ) ):
        thrustingBodyStateFunction_( thrustingBodyStateFunction ),
        centralBodyStateFunction_( centralBodyStateFunction ),
        centralBodyGravitationalParameterFunction_( centralBodyGravitationalParameterFunction ),
        costateFunction_( costateFunction ),
        maximumThrustMagnitude_( thrustMagnitude ),
        specificImpulseFunction_( specificImpulseFunction ),
        thrustingBodyMassFunction_( thrustingBodyMassFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN ),
        customThrustResetFunction_( customThrustResetFunction ){ }

    //! Destructor.
    ~MeeCostatesBangBangThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            currentSpecificImpulse_ = specificImpulseFunction_( time );
            Eigen::VectorXd costates_ = costateFunction_( time );

            Eigen::Vector6d thrustBodyState = thrustingBodyStateFunction_( );
            Eigen::Vector6d centralBodyState = centralBodyStateFunction_( );

            // Get the current state in cartesian coordinates and keplerian elements, and some convenient parameters
            Eigen::Vector6d currentState = thrustBodyState - centralBodyState;
            double centralBodyGravitationalParameter = centralBodyGravitationalParameterFunction_( );

            // Obtain ModifiedEquinoctial elements, flag of 0 indicates that singularity occurs at 180 deg inclination.
            Eigen::Vector6d modifiedEquinoctialElements =
                    orbital_element_conversions::convertCartesianToModifiedEquinoctialElements(
                        currentState, centralBodyGravitationalParameter, 0 );

            double currentMass = thrustingBodyMassFunction_();
            double currentSpecificImpulse =  specificImpulseFunction_( time );

            // Retrieve modified equinoctial elements.
            double p = modifiedEquinoctialElements[ orbital_element_conversions::semiParameterIndex ];
            double f = modifiedEquinoctialElements[ orbital_element_conversions::fElementIndex ];
            double g = modifiedEquinoctialElements[ orbital_element_conversions::gElementIndex ];
            double h = modifiedEquinoctialElements[ orbital_element_conversions::hElementIndex ];
            double k = modifiedEquinoctialElements[ orbital_element_conversions::kElementIndex ];
            double L = modifiedEquinoctialElements[ orbital_element_conversions::trueLongitudeIndex ];

            double w1 = 1.0 + f * std::cos( L ) + g * std::sin( L );
            double w2 = 1.0 + h * h + k * k;

            double lp = costates_[0];
            double lf = costates_[1];
            double lg = costates_[2];
            double lh = costates_[3];
            double lk = costates_[4];
            double lm = costates_[5];

            // Compute all required auxiliary variables to compute optimal angle alpha.
            double lambdaAp = lp * ( 2.0 * p ) / w1;
            double lambdaAf1 = lf * std::sin( L );
            double lambdaAg1 = lg * std::cos( L );
            double lambdaAf2 = lf / w1 * ( ( w1 + 1.0 ) * std::cos( L ) + f );
            double lambdaAg2 = lg / w1 * ( ( w1 + 1.0 ) * std::sin( L ) + g );

            double optimalAlpha = std::atan2(- ( lambdaAf1 - lambdaAg1 ), - ( lambdaAp + lambdaAf2 + lambdaAg2 ));
            double sinOptimalAlpha = sin(optimalAlpha);
            double cosOptimalAlpha = cos(optimalAlpha);

            // Compute all required auxiliary variables to compute optimal angle beta.
            double lambdaBp  = lp * ( 2.0 * p ) / w1 * cosOptimalAlpha;
            double lambdaBf1 = lf * std::sin( L ) * sinOptimalAlpha;
            double lambdaBg1 = lg * std::cos( L ) * sinOptimalAlpha;
            double lambdaBf2 = lf * ( ( 1.0 + w1 ) * std::cos( L ) + f ) / w1 * cosOptimalAlpha;
            double lambdaBg2 = lg * ( ( 1.0 + w1 ) * std::sin( L ) + g ) / w1 * cosOptimalAlpha;
            double lambdaBf3 = lf * ( g / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
            double lambdaBg3 = lg * ( f / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
            double lambdaBh  = lh * ( w2 * std::cos( L ) ) / ( 2.0 * w1 );
            double lambdaBk  = lk * ( w2 * std::sin( L ) ) / ( 2.0 * w1 );

            double optimalBeta = std::atan2(- ( - lambdaBf3 + lambdaBg3 + lambdaBh + lambdaBk ), - ( lambdaBp + lambdaBf1 - lambdaBg1 + lambdaBf2 + lambdaBg2 ));
            double sinOptimalBeta = sin(optimalBeta);
            double cosOptimalBeta = cos(optimalBeta);

            double thrustMagLeft = (1.0 / currentMass) *
                                   (lambdaBp * cosOptimalBeta +
                                    lambdaBh * sinOptimalBeta +
                                    lambdaBk * sinOptimalBeta +
                                    lambdaBf1 * cosOptimalBeta +
                                    lambdaBf2 * cosOptimalBeta -
                                    lambdaBf3 * sinOptimalBeta -
                                    lambdaBg1 * cosOptimalBeta +
                                    lambdaBg2 * cosOptimalBeta +
                                    lambdaBg3 * sinOptimalBeta);

            double thrustMagRight = lm * 1.0 / ( currentSpecificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION);

            // Switching function for the thrust magnitude.
            double thrustMagnitudeSwitchingCondition = thrustMagLeft - thrustMagRight;

            // Compute current thrust magnitude and specific impulse.
            // TODO: Find a better way to handle NaN issues
            if ( isnan(thrustMagnitudeSwitchingCondition) || thrustMagnitudeSwitchingCondition > 0.0 )
            {
                currentThrustMagnitude_ = 0.0;
            }
            else
            {
                currentThrustMagnitude_ = maximumThrustMagnitude_;

            }

            currentTime_ = time;
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustMagnitude( )
    {
        return currentThrustMagnitude_;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( )
    {
        if( currentThrustMagnitude_ != 0.0 )
        {
            return propulsion::computePropellantMassRateFromSpecificImpulse(
                        currentThrustMagnitude_, currentSpecificImpulse_ );
        }
        else
        {
            return 0.0;
        }
    }


    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    {
        if( !( customThrustResetFunction_ == nullptr ) )
        {
            customThrustResetFunction_( currentTime );\
        }
    }

private:

    //!  Function returning the state of the body under thrust as a function of time.
    std::function< Eigen::Vector6d( ) > thrustingBodyStateFunction_;

    //!  Function returning the state of the central body as a function of time.
    std::function< Eigen::Vector6d( ) > centralBodyStateFunction_;

    //!  Function returning the gravitational parameter of the central body as a function of time.
    std::function< double( ) > centralBodyGravitationalParameterFunction_;

    //! General function which gives the costates vector as a function of time.
    std::function< Eigen::VectorXd( const double ) > costateFunction_;

    //! Thrust magnitude when the engine is on.
    double maximumThrustMagnitude_;

    //! Function returning specific impulse as a function of time.
    std::function< double( const double ) > specificImpulseFunction_;

    //! Function returning the mass of the body under thrust as a function of time.
    std::function< double( ) > thrustingBodyMassFunction_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

    std::function< void( const double ) > customThrustResetFunction_;

};


//! Class to compute the engine thrust and mass rate from EngineModel object(s).
class ThrustMagnitudeFromEngineWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor for single engine used for thrust.
    /*!
     * Constructor for single engine used for thrust
     * \param engineModel Engine model used to compute thrust/mass rate.
     */
    ThrustMagnitudeFromEngineWrapper(
            const std::shared_ptr< system_models::EngineModel > engineModel ):
        currentThrust_( TUDAT_NAN ), currentMassRate_( TUDAT_NAN )
    {
        engineModels_.push_back( engineModel );
    }

    //! Constructor for multiple engines used for thrust.
    /*!
     * Constructor for multiple engines used for thrust
     * \param engineModels List of engine models used to compute thrust/mass rate.
     */
    ThrustMagnitudeFromEngineWrapper(
            const std::vector< std::shared_ptr< system_models::EngineModel > > engineModels ):
        engineModels_( engineModels ), currentThrust_( TUDAT_NAN ), currentMassRate_( TUDAT_NAN )
    { }

    //! Destructor.
    ~ThrustMagnitudeFromEngineWrapper( ){ }


    //! Function to update the thrust magnitude to the current time from the engine models.
    /*!
     *  Function to update the thrust magnitude to the current time from the engine models
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ = time ) )
        {
            currentThrust_ = 0.0;
            currentMassRate_ = 0.0;

            // Update engine models
            for( unsigned int i = 0; i < engineModels_.size( ); i++ )
            {
                engineModels_.at( i )->updateEngineModel( time );
            }

            // Add thrust and mass rate.
            for( unsigned int i = 0; i < engineModels_.size( ); i++ )
            {
                currentThrust_ += engineModels_.at( i )->getCurrentThrust( );
                currentMassRate_ += engineModels_.at( i )->getCurrentMassRate( );

            }
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustMagnitude( )
    {
        return currentThrust_;
    }

    //! Function to return the current mass rate
    /*!
     * Function to return the current mass rate, as computed by last call to update member function.
     * \return Current mass rate
     */
    double getCurrentMassRate( )
    {
        return currentMassRate_;
    }

protected:

    //! List of engine models used to compute thrust/mass rate.
    std::vector< std::shared_ptr< system_models::EngineModel > > engineModels_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrust_;

    //! Current mass rate, as computed by last call to update member function.
    double currentMassRate_;
};

//! Variables on which parameterized thrust can depend.
enum ThrustIndependentVariables
{
    mach_number_dependent_thrust,
    altitude_dependent_thrust,
    density_dependent_thrust,
    dynamic_pressure_dependent_thrust,
    pressure_dependent_thrust,
    guidance_input_dependent_thrust,
    throttle_dependent_thrust
};

//! Class to compute the engine thrust and specific impulse as a parameterized function of any number of indepedent
//! variables.
/*!
 *  Class to compute the engine thrust and specific impulse as a parameterized function of any number of indepedent
 *  variables.  The physical meaning of the variables must be defined here, selecting from the options in the
 *  ThrustIndependentVariables enum, and they are automatically retrieved from the relevant environment models during the
 *  propagation.
 *  Note that any number of user-specific functions may be included, as a  guidance_input_dependent_thrust type or
 *  throttle_dependent_thrust. Note that a throttle_dependent_thrust may not be used as one of the independent variables
 *  of the thrust magnitude for this class. This setting is parsed through the ParameterizedThrustMagnitudeSettings
 *  class, which is the class that is typically used to create this class.
 */
class ParameterizedThrustMagnitudeWrapper: public ThrustMagnitudeWrapper
{
public:

    //! Constructor for parameterized thrust and specific impulse.
    /*!
     * Constructor, defines the functions for thrust and specific impulse, as well as the physical meaning of each of the
     * independent variables.
     * \param thrustMagnitudeFunction Function returning the current thrust as a function of the independent variables.
     * \param specificImpulseFunction  Function returning the current specific impulse as a function of the
     * independent variables.
     * \param thrustInputVariableFunctions List of functions returning input variables for the thrust
     * The order of the functions in this vector is passed to the thrustMagnitudeFunction in the same order as
     * entries of this vector.
     * \param specificImpulseInputVariableFunctions List of functions returning input variables for the
     * specific impulse. The order of the functions in this vector is passed to the specificImpulseFunction in the same order
     * as entries of this vector.
     * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
     * the thrustMagnitudeFunction function.
     * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
     * input to the specificImpulseDependentVariables function.
     * \param inputUpdateFunction Function that is called to update the user-defined guidance to the current time
     * (empty by default).
     */
    ParameterizedThrustMagnitudeWrapper(
            const std::function< double( const std::vector< double >& ) > thrustMagnitudeFunction,
            const std::function< double( const std::vector< double >& ) > specificImpulseFunction,
            const std::vector< std::function< double( ) > > thrustInputVariableFunctions,
            const std::vector< std::function< double( ) > > specificImpulseInputVariableFunctions,
            const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
            const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables,
            const std::function< void( const double) > inputUpdateFunction =
            std::function< void( const double) >( ) ):
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        thrustInputVariableFunctions_( thrustInputVariableFunctions ),
        specificImpulseInputVariableFunctions_( specificImpulseInputVariableFunctions ),
        thrustIndependentVariables_( thrustIndependentVariables ),
        specificImpulseDependentVariables_( specificImpulseDependentVariables ),
        inputUpdateFunction_( inputUpdateFunction ),
        currentThrustMagnitude_( TUDAT_NAN ),
        currentSpecificImpulse_( TUDAT_NAN )
    {
        if( thrustInputVariableFunctions_.size( ) != thrustIndependentVariables_.size( ) )
        {
            throw std::runtime_error( "Error in parameterized thrust, inconsistent number of user-defined input variables for thrust" );
        }

        if( specificImpulseInputVariableFunctions_.size( ) != specificImpulseDependentVariables_.size( ) )
        {
            throw std::runtime_error( "Error in parameterized thrust, inconsistent number of user-defined input variables for Isp" );
        }

        currentThrustInputVariables_.resize( thrustInputVariableFunctions.size( ) );
        currentSpecificImpulseInputVariables_.resize( specificImpulseInputVariableFunctions.size( ) );
    }

    //! Destructor.
    ~ParameterizedThrustMagnitudeWrapper( ){ }

    //! Function to update the thrust magnitude to the current time.
    /*!
     *  Function to update the thrust magnitude to the current time.
     *  \param time Time to which the model is to be updated.
     */
    void update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            if( !( inputUpdateFunction_ == nullptr ) )
            {
                inputUpdateFunction_( time );
            }
            // Retrieve thrust independent variables
            for( unsigned int i = 0; i < thrustInputVariableFunctions_.size( ); i++ )
            {
                currentThrustInputVariables_[ i ] = thrustInputVariableFunctions_.at( i )( );
            }

            // Compute thrust
            currentThrustMagnitude_ = thrustMagnitudeFunction_( currentThrustInputVariables_ );

            // Retrieve specific impulse independent variables
            for( unsigned int i = 0; i < specificImpulseInputVariableFunctions_.size( ); i++ )
            {
                currentSpecificImpulseInputVariables_[ i ] = specificImpulseInputVariableFunctions_.at( i )( );
            }

            // Compute specific impulse
            currentSpecificImpulse_ = specificImpulseFunction_( currentSpecificImpulseInputVariables_ );
        }
    }

    //! Function to reset the current time of the thrust model derived class.
    /*!
     *  Function to reset the current time of the thrust model derived class. Function is typically used to reset the time
     *  to NaN, signalling the need for a recomputation of all required quantities.
     *  \param currentTime New current time to be set in model.
     */
    virtual void resetDerivedClassCurrentTime( const double currentTime = TUDAT_NAN )
    {
        if( !( inputUpdateFunction_ == nullptr ) )
        {
            inputUpdateFunction_( currentTime );
        }
    }

    //! Function to return the current thrust magnitude
    /*!
     * Function to return the current thrust magnitude, as computed by last call to update member function.
     * \return Current thrust magnitude
     */
    double getCurrentThrustMagnitude( )
    {
        return currentThrustMagnitude_;
    }

    //! Function to return the current mass rate.
    /*!
     * Function to return the current mass rate, computed from quantities set by last call to update member function.
     * \return Current mass rate.
     */
    double getCurrentMassRate( )
    {
        return propulsion::computePropellantMassRateFromSpecificImpulse(
                    currentThrustMagnitude_, currentSpecificImpulse_ );
    }

private:

    //! Function returning the current thrust as a function of its independent variables
    std::function< double( const std::vector< double >& ) > thrustMagnitudeFunction_;

    //! Function returning the current specific impulse as a function of its independent variables
    std::function< double( const std::vector< double >& ) > specificImpulseFunction_;

    //!  List of functions returning input variables for thrustMagnitudeFunction_.
    std::vector< std::function< double( ) > > thrustInputVariableFunctions_;

    //!  List of functions returning input variables for specificImpulseFunction_.
    std::vector< std::function< double( ) > > specificImpulseInputVariableFunctions_;

    //! List of identifiers for the physical meaning of each of the entries of the input to the thrustMagnitudeFunction
    //! function.
    std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables_;

    //! List of identifiers for the physical meaning of each of the entries of the input to the
    //! specificImpulseDependentVariables function.
    std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables_;

    //! List of current input data to thrust function
    std::vector< double > currentThrustInputVariables_;

    //! List of current input data to specific impulse function.
    std::vector< double > currentSpecificImpulseInputVariables_;


    //! Function that is called to update the user-defined guidance to the current time
    const std::function< void( const double) > inputUpdateFunction_;

    //! Current thrust magnitude, as computed by last call to update member function.
    double currentThrustMagnitude_;

    //! Current specific impulse, as computed by last call to update member function.
    double currentSpecificImpulse_;

};

} // namespace propulsion

} // namespace tudat

#endif // TUDAT_THRUSTMAGNITUDEWRAPPER_H
