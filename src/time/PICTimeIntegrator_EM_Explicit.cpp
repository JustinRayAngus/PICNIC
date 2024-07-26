#include "PICTimeIntegrator_EM_Explicit.H"
#include "System.H"

#include "NamespaceHeader.H"

void PICTimeIntegrator_EM_Explicit::initialize( const std::string&  a_restart_file_name )
{
    m_fields->copyEToVec( m_E );
    m_fields->copyBToVec( m_B );
    m_fields->copyEoldToVec( m_Eold );
    m_fields->copyBoldToVec( m_Bold );
}

int PICTimeIntegrator_EM_Explicit::prepForCheckpoint() const
{
    m_fields->copyEoldFromVec( m_Eold );
    m_fields->copyBoldFromVec( m_Bold );
    return 1;
}

void PICTimeIntegrator_EM_Explicit::preTimeStep(  const Real a_time,
                                                  const Real a_dt,
                                                  const int )
{  
    CH_TIME("PICTimeIntegrator_EM_Explicit::preTimeStep()");
     
    const Real halfDt = 0.5*a_dt;
    const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
    const int num_species = pic_species_ptr_vect.size();

    // initial advance of particles positions by 1/2 time step
    if (a_time==0.0 && m_init_half_advance) {
        for (int sp=0; sp<num_species; sp++) {
            auto species(pic_species_ptr_vect[sp]);
            species->createInflowParticles( a_time, halfDt );
            species->advancePositionsExplicit( halfDt );
        }
    }
  
    // update old particle values and create inflow particles  
    for (int sp=0; sp<num_species; sp++) {
        auto species(pic_species_ptr_vect[sp]);
        species->updateOldParticlePositions();
        species->updateOldParticleVelocities();
        species->createInflowParticles( a_time, a_dt );
    }
  
    if (a_time == 0.0 && m_init_half_advance) {
        m_fields->computeRHSElectricField( halfDt );
        m_fields->copyERHSToVec( m_FE );
        m_E = m_Eold + m_FE;
        m_fields->copyEFromVec( m_E );
        if( m_fields->usePoisson() ) {
           m_pic_species->setChargeDensityOnNodes( m_fields->useFiltering() );
           const LevelData<NodeFArrayBox>& pic_rho = m_pic_species->getChargeDensityOnNodes();
           m_fields->enforceGaussLaw( pic_rho, a_time+0.5*a_dt );
        }
        m_fields->applyBCs_electricField( a_time+0.5*a_dt );
        m_fields->copyEToVec( m_E );
    }

    // update old field values 
    m_Eold = m_E;
    m_Bold = m_B;
    m_fields->setEoldFromE(); // needed for energy diagnostic
    //m_fields->setBoldFromB(); // not needed for energy diagnostic

}

void PICTimeIntegrator_EM_Explicit::timeStep( const Real a_time,
                                              const Real a_dt,
                                              const int )
{  
    CH_TIME("PICTimeIntegrator_EM_Explicit::timeStep()");
  
    // explicit leap-frog time advance of particles and fields
    // Start: we have B_g^n, v_p^n, E_g^{n+1/2}, and x_p^{n+1/2}, time_n
    // End: we have B_g^{n+1}, v_p^{n+1}, E_g^{n+3/2}, and x_p^{n+3/2} 
   
    const Real halfDt = 0.5*a_dt;
    const Real half_time = a_time + halfDt;
    const PicSpeciesPtrVect& pic_species_ptr_vect = m_pic_species->getPtrVect();
    const int num_species = pic_species_ptr_vect.size();
 
    // Step 1: advance B from t_{n} to t_{n+1/2} using E_{n+1/2}
    m_fields->computeRHSMagneticField( halfDt );
    m_fields->copyBRHSToVec( m_FB );
    m_B = m_Bold + m_FB;
    m_fields->updatePhysicalState( m_B, half_time, b_only );
    m_fields->updateBoundaryProbes( a_dt );

    // Step 2: compute Ep and Bp at t_{n+1/2} and xp_{n+1/2}
    //         advance vp from t_{n} to t_{n+1}
    if(m_fields->useFiltering()) { m_fields->setFilteredFields(); }
    for (int sp=0; sp<num_species; sp++) {
        auto species(pic_species_ptr_vect[sp]);
        species->interpolateFieldsToParticles( *m_fields );
        species->addExternalFieldsToParticles( *m_fields );
        if (m_strang_splitting) {
            species->advanceVelocities( halfDt, false );
            species->updateOldParticleVelocities();
        }
        else if (m_average_v_deposit) {
            species->advanceVelocities( a_dt, false );
            species->updateOldParticleVelocities();
            // uncomment lines below to bring xp to t_{n+1/2}
            // prior to scattering
            //species->advancePositionsExplicit( halfDt );
            //species->applyInertialForces( a_dt, false, true, true );
            //species->applyBCs( true, a_time );
        }
        else { species->advanceVelocities( a_dt, false ); }
    }
  
    if (m_post_push_scatter) { m_system->scatterParticles( a_dt ); }

    // advance positions from t_{n+1/2} to t_{n+1} using vp_{n+1}
    const Real time_np1 = a_time + a_dt;
    for (int sp=0; sp<num_species; sp++) {
        auto species(pic_species_ptr_vect[sp]);
        if (m_strang_splitting) {
            species->updateOldParticleVelocities();
            species->advanceVelocities( halfDt, false );
        }
        else if (m_average_v_deposit) {
            species->averageVelocities();
        }
        species->advancePositionsExplicit( halfDt );
        species->applyInertialForces( a_dt, false, true, true );
        species->applyBCs( true, time_np1 );
    }

    // complete advance of B from t_{n+1/2} to t_{n+1}
    m_B = 2.0*m_B - m_Bold;
    m_fields->updatePhysicalState( m_B, time_np1, b_only );

    // compute current density at t_{n+1} and 
    // advance E from t_{n+1/2} to t_{n+3/2} using B_{n+1} and J_{n+1}
    if (m_fields->advanceE()) {
        m_pic_species->setCurrentDensity( *m_fields, a_dt, true, true );
        if (m_fields->useFiltering()) { m_pic_species->filterJ( *m_fields, time_np1 ); }
        const LevelData<EdgeDataBox>& pic_J = m_pic_species->getCurrentDensity();
        const LevelData<NodeFArrayBox>& pic_Jv = m_pic_species->getVirtualCurrentDensity();
        m_fields->setCurrentDensity( pic_J, pic_Jv );
    }
#if CH_SPACEDIM==1
    m_fields->applyAbsorbingBCs( time_np1, a_dt );
#endif
    m_fields->computeRHSElectricField( a_dt );
    m_fields->copyERHSToVec( m_FE );
    m_E = m_Eold + m_FE;
    const Real time_E = time_np1 + halfDt;

    m_fields->copyEFromVec( m_E );
    if(m_fields->usePoisson()) {
        m_pic_species->setChargeDensityOnNodes( m_fields->useFiltering() );
        const LevelData<NodeFArrayBox>& pic_rho = m_pic_species->getChargeDensityOnNodes();
        m_fields->enforceGaussLaw( pic_rho, time_E );
    }
    m_fields->applyBCs_electricField( time_E );
    m_fields->copyEToVec( m_E );

    if (!m_post_push_scatter) { m_system->scatterParticles( a_dt ); }

    // convert positions from t_{n+1} to t_{n+3/2}
    for (int sp=0; sp<num_species; sp++) {
        auto species(pic_species_ptr_vect[sp]);
        species->advancePositions_2ndHalf();
        species->applyBCs( false, time_E );
        if (m_average_v_deposit) { // convert vp from averaged to post-scatter
            species->advanceVelocities_2ndHalf();
        }
    }

}

#include "NamespaceFooter.H"

