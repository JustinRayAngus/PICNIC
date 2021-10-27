#ifndef _PIC_TIME_INTEGRATOR_H_
#define _PIC_TIME_INTEGRATOR_H_

#include <string>
#include "REAL.H"

#include "PicSpecies.H"
#include "ElectroMagneticFields.H"

#include "NamespaceHeader.H"

class System;
typedef std::string PICTIType;

class PICTimeIntegrator
{
  public:

    PICTimeIntegrator() {}

    virtual ~PICTimeIntegrator() {}

    virtual void define( System* const,
                         const PicSpeciesPtrVect&, 
                         ElectroMagneticFields* const ) = 0;

    virtual void initialize() { }

    virtual void preTimeStep( const Real, const Real, const int ) { }

    virtual void timeStep( const Real, const Real, const int ) = 0;

    virtual void postTimeStep( const Real, const Real ) { }

    virtual PICTIType type() const = 0;

    virtual bool isDefined() const = 0;

    virtual void setSolverParams( const Real, const Real, const int ) { }

  protected:

  private:

};

#include "NamespaceFooter.H"
#endif
