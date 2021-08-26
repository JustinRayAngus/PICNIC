#ifndef  _SYSTEM_H_
#define  _SYSTEM_H_

#include <string>
#include <vector>
#include <time.h>

#include "ParmParse.H"
#include "CH_HDF5.H"

#include "PicSpecies.H"
#include "dataFileIO.H"

#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "DomainGrid.H"
#include "CodeUnits.H"
#include "SpecialOperator.H"
#include "ElectroMagneticFields.H"
#include "Scattering.H"

#include "MeshInterp.H"

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

      virtual void initialize( const int cur_step, const double cur_time );
      
      void writePlotFile(const int cur_step, const double cur_time);
      
      void writeHistFile(const int cur_step, const double cur_time, const bool startup_flag);

      Real fieldsDt( const int a_step_number );
      
      Real partsDt( const int a_step_number );
      
      Real scatterDt( const int a_step_number );
      
      Real specialOpsDt( const int a_step_number );
      
      void adaptDt( bool&  a_adapt_dt );

      void preTimeStep( const Real&  a_cur_time, 
                        const Real&  a_dt,
                        const int&   a_step_number );
      
      void advance( Real&  a_cur_time, Real&  a_dt, int&   a_step_number );

      void advance_DSMC( const Real&  a_cur_time, Real&  a_dt );
      void advance_PICMC_EXPLICIT( const Real&  a_cur_time, Real&  a_dt );
      void advance_PICMC_SEMI_IMPLICIT( const Real&  a_cur_time, Real&  a_dt );
      void advance_PICMC_FULLY_IMPLICIT( const Real&  a_cur_time, Real&  a_dt );

      void setChargeDensity();
      void setCurrentDensity();

      void scatterParticles( const Real&  a_dt );
      
      void postTimeStep( const Real&  a_cur_time, 
                         const Real&  a_dt,
                         const int&   a_step_number );
      
   private:

      void parseParameters( ParmParse& ppsys );
      
      void createProblemDomain();
      
      void getDisjointBoxLayout( DisjointBoxLayout& grids );
      
      void setupHistFile();
      
      void createMeshInterp();

      void createState( ParmParse& ppsys );
   
      void createPICspecies();
      
      void createEMfields();
      
      void createScattering();
      
      void createSpecialOperators();
      
      enum ADVANCE_METHOD {
         DSMC,
         PICMC_EXPLICIT,
         PICMC_SEMI_IMPLICIT,
         PICMC_FULLY_IMPLICIT
      };
      int m_iter_max;
      Real m_theta;
     
      ADVANCE_METHOD m_advance_method;

      DomainGrid* m_mesh;
      CodeUnits*  m_units;
      dataFileIO* m_dataFile;
      MeshInterp* m_meshInterp;
     
      PicSpeciesPtrVect m_pic_species_ptr_vect;
      
      ScatteringPtrVect m_scattering_ptr_vect;
     
      RefCountedPtr<ElectroMagneticFields> m_electromagneticFields;
      
      RefCountedPtr<SpecialOperator> m_specialOps;

      bool m_writeSpeciesChargeDensity = false;
      bool m_writeSpeciesCurrentDensity = false;
       
      bool m_verbosity = true;
      bool m_use_scattering = false;
      bool m_use_specialOps = false;

      // problem domain 
      //
      ProblemDomain     m_domain;
      std::string       m_geom_type; 
      int               m_num_ghosts;
      std::vector<int>  m_num_cells;
      std::vector<bool> m_is_periodic;
      std::vector<int>  m_config_decomp;

      Real m_rtol, m_atol;
};

#include "NamespaceFooter.H"

#endif
