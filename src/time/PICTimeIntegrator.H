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

    virtual void initialize( const std::string&  a_restart_file_name ) { }
    
    virtual int prepForCheckpoint() const { return 0; }

    virtual void preTimeStep( const Real, const Real, const int ) { }

    virtual void timeStep( const Real, const Real, const int ) = 0;

    virtual void postTimeStep( const Real, const Real ) { }

    virtual PICTIType type() const = 0;

    virtual bool isDefined() const = 0;

    virtual void printParams() const = 0;
    
    virtual void getConvergenceParams( int&   a_l_exit_status,
                                       int&   a_l_last_iter, 
                                       int&   a_l_total_iter, 
                                       int&   a_nl_exit_status, 
                                       int&   a_nl_iter, 
                                       int&   a_nl_total_iter, 
                                       Real&  a_nl_abs_res, 
                                       Real&  a_nl_rel_res ) const 
    {
      a_l_exit_status = 0;
      a_l_last_iter = 0;
      a_l_total_iter = 0; 
      a_nl_exit_status = 0;
      a_nl_iter = 0;
      a_nl_total_iter = 0;
      a_nl_abs_res = 0.0; 
      a_nl_rel_res = 0.0; 
    }

  protected:

  private:

};

#include "NamespaceFooter.H"
#endif
