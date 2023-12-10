#ifndef _PIC_TIME_INTEGRATOR_EM_THETAIMPLICIT_H_
#define _PIC_TIME_INTEGRATOR_EM_THETAIMPLICIT_H_

#include "PICTimeIntegrator.H"
#include "NonlinearSolverLibrary.H"
#include "ODEVector.H"
#include "EMResidualFunction.H"

#define PIC_EM_THETA_IMPLICIT "PIC_EM_THETA_IMPLICIT"
#define PICMC_FULLY_IMPLICIT "PICMC_FULLY_IMPLICIT"

#include "NamespaceHeader.H"

class PICTimeIntegrator_EM_ThetaImplicit : public PICTimeIntegrator
{
  public:

    PICTimeIntegrator_EM_ThetaImplicit() 
    {
      m_pic_species = nullptr;
      m_fields = nullptr;

      m_nlsolver_type = _NLSOLVER_PICARD_;
      m_nlsolver = nullptr;
      m_func = nullptr;

      m_theta = 0.5;
      m_pc_update_freq = INT_MAX;
      m_pc_update_newton = true;
      m_is_defined = false;
      m_ec_semi_implicit = false;
    }

    virtual ~PICTimeIntegrator_EM_ThetaImplicit() 
    {
      delete m_nlsolver;
      if (m_func) delete m_func;
    }

    virtual 
    void define( PicSpeciesInterface* const,
                 EMFields* const );
    
    virtual void printParams() const;

    virtual void initialize( const std::string&  a_restart_file_name );
    
    inline virtual 
    bool isDefined() const { return m_is_defined; }

    inline virtual 
    PICTIType type() const { return PIC_EM_THETA_IMPLICIT; }

    inline void preTimeStep( const Real  a_time, 
                             const Real  a_dt, 
                             const int   a_step )
    {
      m_Uold = m_U;
      m_nl_iter = 0;

      // update solver verbosity 
      bool this_verbosity = false;
      int cout_step_interval = 1;
      if(a_step<10) cout_step_interval = 1;
      else if(a_step<100) cout_step_interval = 10;
      else cout_step_interval = 100;
      
      if( ((a_step+1) % cout_step_interval)==0 ) {
         this_verbosity = true;
      }
      m_nlsolver->verbose(this_verbosity);
  
      // update old particle values and create inflow particles  
      const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
      const int num_species = pic_species_ptr_vect.size();
      for (int sp=0; sp<num_species; sp++) {
        auto species(pic_species_ptr_vect[sp]);
        species->updateOldParticlePositions();
        species->updateOldParticleVelocities();
        species->createInflowParticles( a_time, a_dt );
        species->injectInflowParticles();
      }

    }

    virtual void timeStep( const Real, const Real, const int );
    
    virtual void preRHSOp( const ODEVector<EMFields>&  a_U,
                           const Real                  a_time, 
                           const Real                  a_dt,
                           const int                   a_nl_iter,
                           const bool                  a_from_emjacobian );
    
    virtual void computeRHS( ODEVector<EMFields>&  a_F,
		       const ODEVector<EMFields>&,
		       const Real                  a_time,
		       const Real                  a_dt,
		       const int                   a_block = -1 );
    
    virtual void updatePhysicalState( ODEVector<EMFields>&  a_U,
                                const int                   a_block,
                                const Real                  a_time );
    
    virtual void updatePrecondMat( BandedMatrix&  a_Pmat,
                             const Real           a_time,
                             const Real           a_dt );
      
    virtual int numPCMatBands() const
    {
      return m_fields->numPCMatBands();
    }

    virtual void getConvergenceParams( int&   a_l_exit_status,
                                       int&   a_l_last_iter, 
                                       int&   a_l_total_iter, 
                                       int&   a_nl_exit_status, 
                                       int&   a_nl_iter, 
                                       int&   a_nl_total_iter, 
                                       Real&  a_nl_abs_res, 
                                       Real&  a_nl_rel_res,
                                       Real&  a_step_wall_time ) const 
    {
       m_nlsolver->getConvergenceParams( a_l_exit_status, a_l_last_iter, a_l_total_iter, 
                                         a_nl_exit_status, a_nl_iter, a_nl_total_iter, 
                                         a_nl_abs_res, a_nl_rel_res );
       a_l_total_iter += m_last_l_total_iter;
       a_nl_total_iter += m_last_nl_total_iter;
       if(a_nl_exit_status==0) { // hack to set probes correct on restart
         a_l_exit_status = m_last_l_exit_status;
         a_l_last_iter = m_last_l_last_iter;
         a_nl_exit_status = m_last_nl_exit_status;
         a_nl_iter = m_last_nl_iter;
         a_nl_abs_res = m_last_nl_abs_res; 
         a_nl_rel_res = m_last_nl_rel_res;
       } 
       a_step_wall_time = m_step_wall_time;
    }

  protected:

    bool m_is_defined;
    bool m_ec_semi_implicit;
    mutable bool m_restart_pc_update_flag = true; // kluge to update pc on restart...
    std::string m_nlsolver_type;


    PicSpeciesInterface*    m_pic_species;
    EMFields*  m_fields;

    Real m_theta;

    ODEVector<EMFields> m_U, m_Uold;
    NonlinearFunction<ODEVector<EMFields>, PICTimeIntegrator>* m_func;
    NonlinearSolver<ODEVector<EMFields>, PICTimeIntegrator>* m_nlsolver;

    int m_pc_update_freq;
    bool m_pc_update_newton;
      
    int m_nl_iter = 0;

  private:

    // below are needed for solver probes on restart
    int m_last_l_exit_status=0; 
    int m_last_nl_exit_status=0;
    int m_last_l_total_iter=0; 
    int m_last_nl_total_iter=0;
    int m_last_l_last_iter=0; 
    int m_last_nl_iter=0;
    Real m_last_nl_abs_res=0.0; 
    Real m_last_nl_rel_res=0.0;
    double m_step_wall_time=0.0;
    double m_pcassembly_walltime=0.0;

};

#include "NamespaceFooter.H"
#endif

