#ifndef  _SYSTEM_H_
#define  _SYSTEM_H_

#include <string>
#include <vector>
#include <time.h>

#include "ParmParse.H"
#include "CH_HDF5.H"

#include "dataFileIO.H"

#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "DomainGrid.H"
#include "CodeUnits.H"

#include "PicSpeciesInterface.H"
#include "ScatteringInterface.H"
#include "ElectroMagneticFields.H"
#include "SpecialOperator.H"

#include "PICTimeIntegratorLibrary.H"

#ifdef CH_USE_TIMER
#include "CH_Timer.H"
#endif

#include "NamespaceHeader.H"

/**
 * System class.
 *
 * This class owns pretty much everything and orchestrates the 
 * user specified time advancement.
 *
 */

class System 
{
   public:

      System( ParmParse&  a_pp );
      virtual ~System();

      virtual void initialize( const int           a_cur_step, 
                               const Real          a_cur_time,
                               const std::string&  a_restart_file_name );
      
      void writePlotFile(const int cur_step, const Real cur_time, const Real cur_dt);
      
      void writeHistFile(const int cur_step, 
                         const Real cur_time, 
                         const Real cur_dt, 
                         const bool startup_flag);
      
      void writeCheckpointFile( HDF5Handle&  a_handle,
                          const int          a_cur_step, 
                          const Real         a_cur_time, 
                          const Real         a_cur_dt );
      
      void readCheckpointFile( const std::string&  a_chkpt_fname );
      
      bool useParts() {
         if (m_pic_species->numSpecies()>0) return true;
         else return false;
      }

      bool useScattering() {
         if (m_scattering->numScatter()>0) return true;
         else return false;
      }
      
      bool useSpecialOps() { return m_use_specialOps; }

      Real fieldsDt( const int a_step_number );
      
      Real partsDt( const int a_step_number );
      
      Real scatterDt( const int a_step_number );
      
      Real specialOpsDt( const int a_step_number );
      
      void adaptDt( bool&  a_adapt_dt );

      void preTimeStep( const Real  a_time, 
                        const Real  a_dt,
                        const int   a_step_number );
      
      void timeStep( const Real, const Real, const int );
      void postTimeStep( Real&, const Real&, int& );

      void setChargeDensity();
      void setCurrentDensity( const bool  a_from_explicit_solver = false );

      void scatterParticles( const Real&  a_dt );
    
      void DefineECSemiImpFlag() { m_ec_semi_implicit = true; }

      inline int nonLinearIter() const { return m_nonlinear_iter; }
      
      void preRHSOp( const ODEVector<System>& a_U,
                     const Real               a_time,
                     const Real               a_dt,
                     const int                a_nlinear_it,
                     const bool               a_from_emjacobian,
                     const EMVecType&         a_type = e_and_b )
      {

        if (a_type == b_only) {
          copyBFromVec( a_U );
          m_emfields->applyBCs_magneticField( a_time );
        } else if (a_type == e_only ) {
          copyEFromVec( a_U );
          m_emfields->applyBCs_electricField( a_time );
        } else {
          copySolutionFromVec( a_U );
          m_emfields->applyBCs_magneticField( a_time );
          m_emfields->applyBCs_electricField( a_time );
        }
        
        if(m_ec_semi_implicit) {
          if(m_nonlinear_iter>1 && !a_from_emjacobian) {
            m_pic_species->postNewtonUpdate( *m_emfields, a_time, a_dt);
          }
          else preRHSOp(a_U, a_dt, a_nlinear_it, a_from_emjacobian);
        }
        else preRHSOp(a_U, a_dt, a_nlinear_it, a_from_emjacobian);

      }
      
      void preRHSOp( const ODEVector<System>&,
                     const Real, 
                     const int,
                     const bool a_from_emjacobian=false);

      void computeRHS(  ODEVector<System>&,
                        const ODEVector<System>&,
                        const Real,
                        const Real,
                        const EMVecType& a_type = e_and_b );

      void computeRHS(  ODEVector<System>&,
                        const ODEVector<System>&,
                        const int,
                        const Real,
                        const Real,
                        const int );

      void updatePhysicalState( ODEVector<System>&,
                          const Real,
                          const EMVecType& a_type = e_and_b );
      
      void updatePhysicalState( ODEVector<System>&,
                          const int,
                          const Real );

      void copyBFromVec( const ODEVector<System>& a_vec, int& );
      void copyEFromVec( const ODEVector<System>& a_vec, int& );
      void copyBToVec( ODEVector<System>& a_vec, int& ) const;
      void copyEToVec( ODEVector<System>& a_vec, int& ) const;
      
      void copyBoldFromVec( const ODEVector<System>& a_vec, int& );
      void copyEoldFromVec( const ODEVector<System>& a_vec, int& );
      void copyBoldToVec( ODEVector<System>& a_vec, int& ) const;
      void copyEoldToVec( ODEVector<System>& a_vec, int& ) const;

      void copyBRHSFromVec( const ODEVector<System>& a_vec, int& );
      void copyERHSFromVec( const ODEVector<System>& a_vec, int& );
      void copyBRHSToVec( ODEVector<System>& a_vec, int& ) const;
      void copyERHSToVec( ODEVector<System>& a_vec, int& ) const;

      void copyBFromVec( const ODEVector<System>& a_vec)
      {
        int offset(0);
        copyBFromVec( a_vec, offset );
      }
      void copyEFromVec( const ODEVector<System>& a_vec)
      {
        int offset(0);
        copyEFromVec( a_vec, offset );
      }
      void copyBToVec( ODEVector<System>& a_vec ) const
      {
        int offset(0);
        copyBToVec( a_vec, offset );
      }
      void copyEToVec( ODEVector<System>& a_vec ) const
      {
        int offset(0);
        copyEToVec( a_vec, offset );
      }

      void copyBoldFromVec( const ODEVector<System>& a_vec)
      {
        int offset(0);
        copyBoldFromVec( a_vec, offset );
      }
      void copyEoldFromVec( const ODEVector<System>& a_vec)
      {
        int offset(0);
        copyEoldFromVec( a_vec, offset );
      }
      void copyBoldToVec( ODEVector<System>& a_vec ) const
      {
        int offset(0);
        copyBoldToVec( a_vec, offset );
      }
      void copyEoldToVec( ODEVector<System>& a_vec ) const
      {
        int offset(0);
        copyEoldToVec( a_vec, offset );
      }

      void copyBRHSFromVec( const ODEVector<System>& a_vec )
      {
        int offset(0);
        copyBRHSFromVec( a_vec, offset );
      }
      void copyERHSFromVec( const ODEVector<System>& a_vec )
      {
        int offset(0);
        copyERHSFromVec( a_vec, offset );
      }
      void copyBRHSToVec( ODEVector<System>& a_vec ) const
      {
        int offset(0);
        copyBRHSToVec( a_vec, offset );
      }
      void copyERHSToVec( ODEVector<System>& a_vec ) const
      {
        int offset(0);
        copyERHSToVec( a_vec, offset );
      }

      void copySolutionFromVec( const ODEVector<System>& a_vec );
      void copySolutionToVec( ODEVector<System>& a_vec) const;
      void copyRHSFromVec( const ODEVector<System>& a_vec);
      void copyRHSToVec( ODEVector<System>& a_vec) const;
      
      inline int numPicardBlocks() const
      {
        return m_emfields->numPicardBlocks();
      }

      inline int getVectorSize() const
      {
        int retval(0);
        if (m_advance_method==PIC_EM_SEMI_IMPLICIT) {
          retval += m_emfields->getVectorSize(e_only);
        } else if (m_advance_method==PIC_EM_THETA_IMPLICIT) {
          retval += m_emfields->getVectorSize(e_and_b);
        }
        return retval;
      }

      inline int getVectorSize( const EMVecType& a_type ) const
      {
        return m_emfields->getVectorSize(a_type);
      }

      inline void updatePrecondMat( BandedMatrix& a_Pmat,
                                    const Real    a_time,
                                    const Real    a_dt )
      {
        CH_TIME("System::updatePrecondMat()");

        a_Pmat.setToIdentityMatrix();
        const Real cnormDt = a_dt*m_units->CvacNorm();
        //use_mm = m_pic_species->useMassMatrices();
        if(m_pic_species->getSigmaxx().isDefined()) {
	   m_pic_species->setMassMatricesForPC(*m_emfields);
	}
        m_emfields->assemblePrecondMatrix( a_Pmat,
                                           m_pic_species->getSigmaxx().isDefined(),
                                           cnormDt );
        a_Pmat.finalAssembly();
        //a_Pmat.writeToFile("pc_matrix.txt",1);
      }

      inline int numPCMatBands() const
      {
        return m_emfields->numPCMatBands();
      }

   private:

      void parseParameters( ParmParse& ppsys );
      
      void createProblemDomain();
      
      void getDisjointBoxLayout( DisjointBoxLayout& grids );
      
      void setupHistFile( const int  a_cur_step );
      
      void createState();
      
      void createEMfields();
      
      void createSpecialOperators();

      PICTIType m_advance_method;
      bool m_implicit_advance;

      DomainGrid* m_mesh;
      CodeUnits*  m_units;
      dataFileIO* m_dataFile;
      PicSpeciesInterface* m_pic_species;
      ScatteringInterface* m_scattering;
     
      RefCountedPtr<ElectroMagneticFields> m_emfields;
      RefCountedPtr<SpecialOperator> m_specialOps;

      bool m_solver_probes;
      bool m_field_probes;
      bool m_scattering_probes;
      std::vector<bool> m_species_probes;
      bool m_field_bdry_probes;
      bool m_species_bdry_probes;
      bool m_species_solver_probes;
      int m_history_precision;      

      bool m_verbosity;
      bool m_use_specialOps;

      int m_nonlinear_iter = 0;
      bool m_ec_semi_implicit = false;

      // problem domain 
      ProblemDomain     m_domain;
      int               m_num_ghosts;
      std::vector<int>  m_num_cells;
      std::vector<bool> m_is_periodic;
      std::vector<int>  m_config_decomp;

      PICTimeIntegrator* m_time_integrator;
};

#include "NamespaceFooter.H"

#endif
