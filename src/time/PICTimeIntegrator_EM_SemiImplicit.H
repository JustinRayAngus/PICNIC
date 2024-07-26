#ifndef _PIC_TIME_INTEGRATOR_EM_SEMI_IMPLICIT_H_
#define _PIC_TIME_INTEGRATOR_EM_SEMI_IMPLICIT_H_

#include "PICTimeIntegrator.H"
#include "NonlinearSolverLibrary.H"
#include "ODEVector.H"

#define PIC_EM_SEMI_IMPLICIT "PIC_EM_SEMI_IMPLICIT"
#define PICMC_SEMI_IMPLICIT "PICMC_SEMI_IMPLICIT"

#include "NamespaceHeader.H"

class PICTimeIntegrator_EM_SemiImplicit : public PICTimeIntegrator
{
  public:

    PICTimeIntegrator_EM_SemiImplicit() 
    {
      m_pic_species = nullptr;
      m_fields = nullptr;
      
      m_nlsolver_type = _NLSOLVER_PICARD_;
      m_nlsolver = nullptr;
      m_func = nullptr;

      m_pc_update_freq = INT_MAX;
      m_pc_update_newton = true;

      m_is_defined = false;
    }

    virtual ~PICTimeIntegrator_EM_SemiImplicit() 
    {
      delete m_nlsolver;
      if (m_func) delete m_func;
    }

    virtual 
    void define( System* const,
                 PicSpeciesInterface* const,
                 EMFields* const );

    virtual void initialize( const std::string&  a_restart_file_name );
    
    virtual void printParams() const;
    
    virtual int prepForCheckpoint() const;

    inline virtual 
    bool isDefined() const { return m_is_defined; }

    inline virtual 
    PICTIType type() const { return PIC_EM_SEMI_IMPLICIT; }

    virtual void preTimeStep( const Real, const Real, const int );
    
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
		       const int );
    
    virtual void updatePhysicalState( ODEVector<EMFields>&  a_U,
                                const int,
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
    std::string m_nlsolver_type;

    int m_pc_update_freq;
    bool m_pc_update_newton;
      
    int m_nl_iter = 0;
    mutable bool m_restart_pc_update_flag = true; // kluge to update pc on restart...


    PicSpeciesInterface*    m_pic_species;
    EMFields*  m_fields;

    ODEVector<EMFields> m_E, m_Eold;
    ODEVector<EMFields> m_B, m_Bold, m_FB;
    NonlinearFunction<ODEVector<EMFields>, PICTimeIntegrator>* m_func;
    NonlinearSolver<ODEVector<EMFields>, PICTimeIntegrator>* m_nlsolver;

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

