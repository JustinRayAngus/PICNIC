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
      m_system = nullptr;
      m_particles.clear();

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
    void define( System* const,
                 const PicSpeciesPtrVect&,
                 ElectroMagneticFields* const );
    
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
      for (int s=0; s<m_particles.size(); s++) {
        auto this_picSpecies(m_particles[s]);
        this_picSpecies->updateOldParticlePositions();
        this_picSpecies->updateOldParticleVelocities();
        this_picSpecies->createInflowParticles( a_time, a_dt );
        this_picSpecies->injectInflowParticles();
      }

    }

    virtual void timeStep( const Real, const Real, const int );
    virtual void postTimeStep( const Real, const Real );

    virtual void getConvergenceParams( int&   a_l_exit_status,
                                       int&   a_l_last_iter, 
                                       int&   a_l_total_iter, 
                                       int&   a_nl_exit_status, 
                                       int&   a_nl_iter, 
                                       int&   a_nl_total_iter, 
                                       Real&  a_nl_abs_res, 
                                       Real&  a_nl_rel_res ) const 
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
    }

  protected:

    bool m_is_defined;
    bool m_ec_semi_implicit;
    mutable bool m_restart_pc_update_flag = true; // kluge to update pc on restart...
    std::string m_nlsolver_type;


    System*                 m_system;
    PicSpeciesPtrVect       m_particles;

    Real m_theta;

    ODEVector<System> m_U, m_Uold;
    NonlinearFunction<ODEVector<System>, System>* m_func;
    NonlinearSolver<ODEVector<System>, System>* m_nlsolver;

    int m_pc_update_freq;
    bool m_pc_update_newton;

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

};

#include "NamespaceFooter.H"
#endif

