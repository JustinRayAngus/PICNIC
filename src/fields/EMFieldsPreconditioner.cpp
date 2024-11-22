#include <map>
#include "EMFields.H"

#include "NamespaceHeader.H"

template<typename kt, typename dt>
static void insertOrAdd( std::map<kt,dt>& a_map,
                         const kt& a_key,
                         const dt& a_value )
{
    auto it = a_map.find(a_key);
    if (it != a_map.end()) {
        a_map[a_key] += a_value;
    } else {
        a_map[a_key] = a_value;
    }

    return;
}

void EMFields::computePrecondMatrixNNZ( const EMVecType& a_vec_type)
{
   /* setting number of non-zero elements per row in preconditioner matrix */

   if (m_pc_diag_only) {
     m_pcmat_nnz = 1;
   } else {
#if CH_SPACEDIM==1
     static int nbands_J_mass_matrix = 1 + 2*3*(m_pc_mass_matrix_width+1);
#elif CH_SPACEDIM==2
     static int nbands_J_mass_matrix =    3
                                        * (1+2*m_pc_mass_matrix_width)
                                        * (1+2*m_pc_mass_matrix_width);
#endif
     if (a_vec_type == e_and_b) {
        static int nbands_EMfields = 1 + 2*SpaceDim;
        m_pcmat_nnz = nbands_EMfields + nbands_J_mass_matrix;
     } else if (a_vec_type == e_only) {
        m_pcmat_nnz = nbands_J_mass_matrix;
     } else if (a_vec_type == curl2) {
        // the following overallocates
        static int nbands_EMfields = 1 + 2*SpaceDim + 2*(SpaceDim-1)*(SpaceDim-1); // too many for 3D
        m_pcmat_nnz = nbands_EMfields + nbands_J_mass_matrix;
     }
   }
}

void EMFields::assemblePrecondMatrix( BandedMatrix&  a_P,
                                const bool           a_use_mass_matrices,
                                const Real           a_dt,
                                const EMVecType&     a_vec_type )
{
  CH_TIME("EMFields::assemblePrecondMatrix()");
  if (a_vec_type == e_only) {
    if (m_pc_diag_only) {
      assemblePrecondMatrixEDiagOnly( a_P, a_use_mass_matrices, a_dt );
    } else {
      assemblePrecondMatrixE( a_P, false, a_use_mass_matrices, a_dt );
    }
  } else if (a_vec_type == e_and_b) {
    if (m_pc_diag_only) {
      assemblePrecondMatrixBDiagOnly( a_P );
      assemblePrecondMatrixEDiagOnly( a_P, a_use_mass_matrices, a_dt );
    } else {
      assemblePrecondMatrixB( a_P, true, a_dt );
      assemblePrecondMatrixE( a_P, true, a_use_mass_matrices, a_dt );
    }
  } else if (a_vec_type == curl2) {
    if (m_pc_diag_only) {
      //TODO
      MayDay::Error("Not yet implemented");
      //assemblePrecondMatrixCurl2DiagOnly( a_P, a_use_mass_matrices, a_dt );
    } else {
      assemblePrecondMatrixCurl2( a_P, true, a_use_mass_matrices, a_dt );
    }
  } else {
      MayDay::Error("EMFields::assemblePrecondMatrix(): not yet implemented");
  }
}

void EMFields::assemblePrecondMatrixBDiagOnly(  BandedMatrix&  a_P )
{
  CH_TIME("EMFields::assemblePrecondMatrixBDiagOnly()");

  const LevelData<FluxBox>& pmap_B_lhs = m_gdofs.dofDataLHSB();
  const LevelData<FArrayBox>& pmap_Bv_lhs = m_gdofs.dofDataLHSBv();

  auto grids(m_mesh.getDBL());
  auto phys_domain( grids.physDomain() );
  const int ghosts(m_mesh.ghosts());

  for (auto dit( grids.dataIterator() ); dit.ok(); ++dit) {

    /* magnetic field */
    for (int dir = 0; dir < SpaceDim; dir++) {
      const FArrayBox& pmap( pmap_B_lhs[dit][dir] );
      auto box( grow(pmap.box(), -ghosts) );
      int idir_bdry = phys_domain.domainBox().bigEnd(dir);
      if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
      for (BoxIterator bit(box); bit.ok(); ++bit) {
        CH_assert(pmap.nComp() == 1);
        int icol = (int) pmap(bit(),0);
        Real val = 1.0;
        a_P.setRowValues( icol, 1, &icol, &val );
      }
    }

#if CH_SPACEDIM<3
    /* virtual magnetic field */
    {
      const FArrayBox& pmap( pmap_Bv_lhs[dit] );
      auto box( grow(pmap.box(), -ghosts) );
      for (BoxIterator bit(box); bit.ok(); ++bit) {
        for (int n(0); n < pmap.nComp(); n++) {
          int icol = (int) pmap(bit(),n);
          Real val = 1.0;
          a_P.setRowValues( icol, 1, &icol, &val );
        }
      }
    }
#endif

  }

  return;
}

void EMFields::assemblePrecondMatrixB(  BandedMatrix&  a_P,
                                  const bool           a_include_EM,
                                  const Real           a_dt )
{
  CH_TIME("EMFields::assemblePrecondMatrixB()");

  auto grids(m_mesh.getDBL());

  // copy over the mass matrix elements used in the PC and do addOp exchange()
  const Real cnormDt = a_dt*m_cvacNorm;

  const LevelData<FluxBox>& pmap_B_lhs = m_gdofs.dofDataLHSB();
  const LevelData<FArrayBox>& pmap_Bv_lhs = m_gdofs.dofDataLHSBv();
#if CH_SPACEDIM>1
  const LevelData<EdgeDataBox>& pmap_E_rhs = m_gdofs.dofDataRHSE();
#endif
  const LevelData<NodeFArrayBox>& pmap_Ev_rhs = m_gdofs.dofDataRHSEv();

#if CH_SPACEDIM==1
  const LevelData<FArrayBox>& Xcc(m_mesh.getXcc());
#endif
  const LevelData<NodeFArrayBox>& Xnc(m_mesh.getXnc());
#if CH_SPACEDIM==2
  const LevelData<FluxBox>& Xfc(m_mesh.getXfc());
#endif

  auto phys_domain( grids.physDomain() );
  const int ghosts(m_mesh.ghosts());
  const RealVect& dX(m_mesh.getdX());

#if CH_SPACEDIM==1
  const string& geom_type = m_mesh.geomType();
#endif

  for (auto dit( grids.dataIterator() ); dit.ok(); ++dit) {

    /* magnetic field */
    for (int dir = 0; dir < SpaceDim; dir++) {

      const FArrayBox& pmap( pmap_B_lhs[dit][dir] );

      auto box( grow(pmap.box(), -ghosts) );
      int idir_bdry = phys_domain.domainBox().bigEnd(dir);
      if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);

      for (BoxIterator bit(box); bit.ok(); ++bit) {

        auto ic = bit();
        CH_assert(pmap.nComp() == 1);
        int pc = (int) pmap(ic, 0);

        int ncols = m_pcmat_nnz;
        std::vector<int> icols(ncols, -1);
        std::vector<Real> vals(ncols, 0.0);
        int ix = 0;

        icols[ix] = pc;
        vals[ix] = 1.0;
        ix++;

#if CH_SPACEDIM==2
        if (a_include_EM) {

          int tdir = (dir + 1) % SpaceDim;
          int sign = 1 - 2*dir;
          if(m_mesh.anticyclic()) sign *= -1;

          IntVect iL(ic);
          Real aL = (-1) * sign * (-cnormDt/dX[tdir]);
          aL *= (m_PC_mask_Ev[dit](iL,2*tdir) * m_advanceB_comp[dir]);
          int pL = (int) pmap_Ev_rhs[dit](iL, 0);

          if (pL >= 0) {
            icols[ix] = pL;
            if (m_mesh.axisymmetric() && (dir==1)) {
              Real local_Xfc = Xfc[dit][dir](ic,0);
              if (local_Xfc==0.0) vals[ix] = 0.0;
              else vals[ix] = -aL*Xnc[dit](iL,0)/local_Xfc;
            } else {
              vals[ix] = -aL;
            }
            ix++;
          }

          IntVect iR(ic); iR[tdir]++;
          Real aR = sign * (-cnormDt/dX[tdir]);
          aR *= (m_PC_mask_Ev[dit](iR,2*tdir+1) * m_advanceB_comp[dir]);
          int pR = (int) pmap_Ev_rhs[dit](iR, 0);

          if (pR >= 0) {
            icols[ix] = pR;
            if (m_mesh.axisymmetric() && (dir==1)) {
              Real local_Xfc = Xfc[dit][dir](ic,0);
              if (local_Xfc==0.0) vals[ix] = 0.0;
              else vals[ix] = -aR*Xnc[dit](iR,0)/local_Xfc;
            } else {
              vals[ix] = -aR;
            }
            ix++;
          }
        }
#endif

        CH_assert(ix <= m_pcmat_nnz);
        CH_assert(ix <= a_P.getNBands());

        a_P.setRowValues( pc, ix, icols.data(), vals.data() );
      }
    }

#if CH_SPACEDIM<3
    /* virtual magnetic field */
    {
      const FArrayBox& pmap( pmap_Bv_lhs[dit] );
      auto box( grow(pmap.box(), -ghosts) );

      for (BoxIterator bit(box); bit.ok(); ++bit) {

        auto ic = bit();

        for (int n(0); n < pmap.nComp(); n++) {
          int pc = (int) pmap(ic, n);

          int ncols = m_pcmat_nnz;
          std::vector<int> icols(ncols, -1);
          std::vector<Real> vals(ncols, 0.0);
          int ix = 0;

          icols[ix] = pc;
          vals[ix] = 1.0;
          ix++;

          if (a_include_EM) {
#if CH_SPACEDIM==1
            int dirj = 0, dirk = 1;
            if (n == 0) {

              int sign = -1;

	          // axisymmetric modifications are needed in 1D sph geom
              // sph: dBth_{i+1/2}/dt = (r_{i+1}*Ephi_{i+1} - r_{i}*Ephi_{i})/dr/r_{i+1/2}

              Real aL = (-1) * sign * (-cnormDt/dX[0]) * m_advanceB_comp[n+1];
              IntVect iL(ic);
              int pL = (int) pmap_Ev_rhs[dit](iL, dirk);

              if (pL >= 0) {
                icols[ix] = pL;
                if(m_mesh.axisymmetric() && geom_type=="sph_R") {
                  vals[ix] = -aL*Xnc[dit](iL,0)/Xcc[dit](ic,0);
                }
                else { vals[ix] = -aL; }
                ix++;
              }

              Real aR =  sign * (-cnormDt/dX[0]) * m_advanceB_comp[n+1];
              IntVect iR(ic); iR[0]++;
              int pR = (int) pmap_Ev_rhs[dit](iR, dirk);

              if (pR >= 0) {
                icols[ix] = pR;
                if(m_mesh.axisymmetric() && geom_type=="sph_R") {
                  vals[ix] = -aR*Xnc[dit](iR,0)/Xcc[dit](ic,0);
                }
                else { vals[ix] = -aR; }
                ix++;
              }

            } else if (n == 1) {

              int sign = 1;

	          // axisymmetric modifications are needed for 1D cyl and 1D sph geoms
              // cyl: dBz_{i+1/2}/dt = -(r_{i+1}*Eth_{i+1} - r_{i}*Eth_{i})/dr/r_{i+1/2}
              // sph: dBphi_{i+1/2}/dt = -(r_{i+1}*Eth_{i+1} - r_{i}*Eth_{i})/dr/r_{i+1/2}

              Real aL = (-1) * sign * (-cnormDt/dX[0]) * m_advanceB_comp[n+1];
              IntVect iL(ic);
              int pL = (int) pmap_Ev_rhs[dit](iL, dirj);

              if (pL >= 0) {
                icols[ix] = pL;
                if(m_mesh.axisymmetric()) {
                  vals[ix] = -aL*Xnc[dit](iL,0)/Xcc[dit](ic,0);
                }
                else { vals[ix] = -aL; }
                ix++;
              }

              Real aR =  sign * (-cnormDt/dX[0]) * m_advanceB_comp[n+1];
              IntVect iR(ic); iR[0]++;
              int pR = (int) pmap_Ev_rhs[dit](iR, dirj);

              if (pR >= 0) {
                icols[ix] = pR;
                if(m_mesh.axisymmetric()) {
                  vals[ix] = -aR*Xnc[dit](iR,0)/Xcc[dit](ic,0);
                }
                else { vals[ix] = -aR; }
                ix++;
              }

            }
#elif CH_SPACEDIM==2
            int dirj = 0, dirk = 1;
            {
              int sign = 1;
              if(m_mesh.anticyclic()) sign *= -1;

              IntVect iL(ic);
              Real aL = (-1) * sign * (-cnormDt/dX[dirj]);
              aL *= (m_PC_mask_E[dit][dirk](iL,0) * m_advanceB_comp[2]);
              int pL = (int) pmap_E_rhs[dit][dirk](iL, 0);

              if (pL >= 0) {
                icols[ix] = pL;
                vals[ix] = -aL;
                ix++;
              }

              IntVect iR(ic); iR[dirj]++;
              Real aR = sign * (-cnormDt/dX[dirj]);
              aR *= (m_PC_mask_E[dit][dirk](iR,1) * m_advanceB_comp[2]);
              int pR = (int) pmap_E_rhs[dit][dirk](iR, 0);

              if (pR >= 0) {
                icols[ix] = pR;
                vals[ix] = -aR;
                ix++;
              }
            }
            {
              int sign = -1;
              if(m_mesh.anticyclic()) sign *= -1;

              IntVect iL(ic);
              Real aL = (-1) * sign * (-cnormDt/dX[dirk]);
              aL *= (m_PC_mask_E[dit][dirj](iL,0) * m_advanceB_comp[2]);
              int pL = (int) pmap_E_rhs[dit][dirj](iL, 0);

              if (pL >= 0) {
                icols[ix] = pL;
                vals[ix] = -aL;
                ix++;
              }

              IntVect iR(ic); iR[dirk]++;
              Real aR = sign * (-cnormDt/dX[dirk]);
              aR *= (m_PC_mask_E[dit][dirj](iR,1) * m_advanceB_comp[2]);
              int pR = (int) pmap_E_rhs[dit][dirj](iR, 0);

              if (pR >= 0) {
                icols[ix] = pR;
                vals[ix] = -aR;
                ix++;
              }
            }
#endif
          }
          CH_assert(ix <= m_pcmat_nnz);
          CH_assert(ix <= a_P.getNBands());

          a_P.setRowValues( pc, ix, icols.data(), vals.data() );
        }
      }
    }
#endif

  }

  return;
}

void EMFields::assemblePrecondMatrixEDiagOnly(  BandedMatrix&  a_P,
                                          const bool           a_use_mass_matrices,
                                          const Real           a_dt )
{
  CH_TIME("EMFields::assemblePrecondMatrixEDiagOnly()");

  auto grids(m_mesh.getDBL());

  // copy over the mass matrix elements used in the PC and do addOp exchange()
  const Real cnormDt = a_dt*m_cvacNorm;
  const Real sig_normC = cnormDt*m_Jnorm_factor;
  int diag_comp_inPlane=0, diag_comp_virtual=0;
  if(a_use_mass_matrices) {
    diag_comp_inPlane = (m_sigma_xx_pc.nComp()-1)/2;
#if CH_SPACEDIM<3
    diag_comp_virtual = (m_sigma_zz_pc.nComp()-1)/2;
#endif
  }

  const LevelData<EdgeDataBox>& pmap_E_lhs = m_gdofs.dofDataLHSE();
  const LevelData<NodeFArrayBox>& pmap_Ev_lhs = m_gdofs.dofDataLHSEv();

#if CH_SPACEDIM<3
  const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();
#endif
  const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();

  auto phys_domain( grids.physDomain() );
  const int ghosts(m_mesh.ghosts());

#define FILTER_PC_TEST

  if (a_use_mass_matrices && m_use_filtering) {
    SpaceUtils::exchangeEdgeDataBox(m_sigma_xx_pc);
#if CH_SPACEDIM==1
    SpaceUtils::exchangeNodeFArrayBox(m_sigma_yy_pc);
#endif
    SpaceUtils::exchangeNodeFArrayBox(m_sigma_zz_pc);
  }

  for (auto dit( grids.dataIterator() ); dit.ok(); ++dit) {

    /* electric field */
    for (int dir = 0; dir < SpaceDim; dir++) {

      const FArrayBox& pmap( pmap_E_lhs[dit][dir] );

      auto box( grow(pmap.box(), -ghosts) );
      for (int adir=0; adir<SpaceDim; ++adir) {
         if (adir != dir) {
            int idir_bdry = phys_domain.domainBox().bigEnd(adir);
            if (box.bigEnd(adir) < idir_bdry) box.growHi(adir, -1);
         }
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {
        auto ic = bit();
        CH_assert( pmap.nComp() == 1);
	    const Real sig_norm_ec = sig_normC/Jec[dit][dir](ic,0);

        int icol = ((int) pmap(ic,0));
        Real val = 1.0;

        if(a_use_mass_matrices) { // dJ_{ic}/dE_{ic}
#ifdef FILTER_PC_TEST
           if(m_use_filtering) {
             val += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-1);
             val += sig_norm_ec*4.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);
             val += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+1);
	       } else {
#endif
             val += sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);
#ifdef FILTER_PC_TEST
	       }
#endif
        }

        a_P.setRowValues( icol, 1, &icol, &val );
      }
    }

#if CH_SPACEDIM<3
    /* virtual electric field */
    {
      const NodeFArrayBox& pmap( pmap_Ev_lhs[dit] );
      auto box( surroundingNodes( grow(pmap.box(), -ghosts) ) );
      for (int dir=0; dir<SpaceDim; ++dir) {
         int idir_bdry = phys_domain.domainBox().bigEnd(dir);
         if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {
        auto ic = bit();
        const Real sig_norm_nc = sig_normC/Jnc[dit](ic,0);
        for (int n(0); n < pmap.nComp(); n++) {
          int icol = (int) pmap(ic, n);
          Real val = 1.0;
#if CH_SPACEDIM==1
          if(a_use_mass_matrices && m_advanceE_comp[1+n]) {
              val += sig_norm_nc * ( n==0 ? m_sigma_yy_pc[dit](ic,diag_comp_virtual)
                                          : m_sigma_zz_pc[dit](ic,diag_comp_virtual) );
          }
#elif CH_SPACEDIM==2
          if(a_use_mass_matrices && m_advanceE_comp[2]) {
            val += sig_norm_nc*m_sigma_zz_pc[dit](ic,diag_comp_virtual);
          }
#endif
          a_P.setRowValues( icol, 1, &icol, &val );
        }
      }
    }
#endif

  }

  return;
}

void EMFields::assemblePrecondMatrixE(  BandedMatrix&  a_P,
                                  const bool           a_include_EM,
                                  const bool           a_use_mass_matrices,
                                  const Real           a_dt )
{
  CH_TIME("EMFields::assemblePrecondMatrixE()");

  auto grids(m_mesh.getDBL());

  // copy over the mass matrix elements used in the PC and do addOp exchange()
  const Real cnormDt = a_dt*m_cvacNorm;
  const Real sig_normC = cnormDt*m_Jnorm_factor;
  int diag_comp_inPlane=0, diag_comp_virtual=0;
  if(a_use_mass_matrices) {
    diag_comp_inPlane = (m_sigma_xx_pc.nComp()-1)/2;
#if CH_SPACEDIM<3
    diag_comp_virtual = (m_sigma_zz_pc.nComp()-1)/2;
#endif
  }

  const LevelData<EdgeDataBox>& pmap_E_lhs = m_gdofs.dofDataLHSE();
  const LevelData<NodeFArrayBox>& pmap_Ev_lhs = m_gdofs.dofDataLHSEv();
  const LevelData<FArrayBox>& pmap_Bv_rhs = m_gdofs.dofDataRHSBv();
  const LevelData<EdgeDataBox>& pmap_E_rhs = m_gdofs.dofDataRHSE();
  const LevelData<NodeFArrayBox>& pmap_Ev_rhs = m_gdofs.dofDataRHSEv();
#if CH_SPACEDIM>1
  const LevelData<FluxBox>& pmap_B_rhs = m_gdofs.dofDataRHSB();
#endif

  const LevelData<FArrayBox>& Xcc(m_mesh.getXcc());
#if CH_SPACEDIM==1
  const LevelData<NodeFArrayBox>& Xnc(m_mesh.getXnc());
#endif
#if CH_SPACEDIM==2
  const LevelData<EdgeDataBox>& Xec(m_mesh.getXec());
#endif

#if CH_SPACEDIM<3
  const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();
#endif
  const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();

  auto phys_domain( grids.physDomain() );
  const int ghosts(m_mesh.ghosts());
  const RealVect& dX(m_mesh.getdX());

#define FILTER_PC_TEST

  if (a_use_mass_matrices && m_use_filtering) {
    SpaceUtils::exchangeEdgeDataBox(m_sigma_xx_pc);
    if(m_pc_mass_matrix_include_ij) {
      SpaceUtils::exchangeEdgeDataBox(m_sigma_xy_pc);
      SpaceUtils::exchangeEdgeDataBox(m_sigma_xz_pc);
    }
#if CH_SPACEDIM==1
    SpaceUtils::exchangeNodeFArrayBox(m_sigma_yy_pc);
    if(m_pc_mass_matrix_include_ij) {
      SpaceUtils::exchangeNodeFArrayBox(m_sigma_yx_pc);
      SpaceUtils::exchangeNodeFArrayBox(m_sigma_yz_pc);
    }
#endif
    if(m_pc_mass_matrix_include_ij) {
      SpaceUtils::exchangeNodeFArrayBox(m_sigma_zx_pc);
      SpaceUtils::exchangeNodeFArrayBox(m_sigma_zy_pc);
    }
    SpaceUtils::exchangeNodeFArrayBox(m_sigma_zz_pc);
  }

#if CH_SPACEDIM==1
  const string& geom_type = m_mesh.geomType();
#endif

  for (auto dit( grids.dataIterator() ); dit.ok(); ++dit) {

    /* electric field */
    for (int dir = 0; dir < SpaceDim; dir++) {

      const FArrayBox& pmap( pmap_E_lhs[dit][dir] );

      auto box( grow(pmap.box(), -ghosts) );
      for (int adir=0; adir<SpaceDim; ++adir) {
         if (adir != dir) {
            int idir_bdry = phys_domain.domainBox().bigEnd(adir);
            if (box.bigEnd(adir) < idir_bdry) box.growHi(adir, -1);
         }
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {

        auto ic = bit();
        CH_assert( pmap.nComp() == 1);
        int pc = (int) pmap(ic, 0);

	      const Real sig_norm_ec = sig_normC/Jec[dit][dir](ic,0);

        int ncols = m_pcmat_nnz;
        std::vector<int> icols(ncols, -1);
        std::vector<Real> vals(ncols, 0.0);

        int ix = 0;
        icols[ix] = pc;
        vals[ix] = 1.0;
        if(a_use_mass_matrices) { // dJ_{ic}/dE_{ic}
#ifdef FILTER_PC_TEST
           if(m_use_filtering) {
             vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-1);
             vals[ix] += sig_norm_ec*4.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);
             vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+1);
             IntVect iL(ic); iL[0] -= 1;
	           int pL = (int) pmap_E_rhs[dit][dir](iL,0);
	           if(pL>=0) {
		           vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iL,diag_comp_inPlane);
		           vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](iL,diag_comp_inPlane+1);
		           vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iL,diag_comp_inPlane+2);
	           }
             IntVect iR(ic); iR[0] += 1;
	           int pR = (int) pmap_E_rhs[dit][dir](iR,0);
	           if(pR>=0) {
		           vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iR,diag_comp_inPlane-2);
		           vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](iR,diag_comp_inPlane-1);
		           vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iR,diag_comp_inPlane);
	           }
	         }
	         else {
             vals[ix] += sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);
	         }
#else
           vals[ix] += sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);
#endif
        }
        ix++;

#if CH_SPACEDIM==1
        if (a_use_mass_matrices && m_advanceE_comp[0]) {

          /* sigma_xx contributions */
          for (int n = 1; n < m_pc_mass_matrix_width+1; n++) {

            if (n > (m_sigma_xx_pc.nComp()-1)/2) break;

            IntVect iL(ic); iL[0] -= n;
            int pL = (int) pmap_E_rhs[dit][dir](iL, 0);
            if (pL >= 0) {
              icols[ix] = pL;
              vals[ix] = 0.0;
#ifdef FILTER_PC_TEST
              if(m_use_filtering && n==1) { // dJ_{ic}/dE_{ic-1}
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iL,diag_comp_inPlane-1);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](iL,diag_comp_inPlane);
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iL,diag_comp_inPlane+1);

                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-2);
                vals[ix] += sig_norm_ec*4.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-1);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);

                IntVect ip(ic); ip[0] += 1;
                int pR = (int) pmap_E_rhs[dit][dir](ip, 0);
                if(pR >= 0) {
                  vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane-2);
                  vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane-1);
                }
              }
              else if(m_use_filtering && n==2) { // dJ_{ic}/dE_{ic-2}
                IntVect im(ic); im[0] -= 1;
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane-2);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane-1);
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane);

                vals[ix] += sig_norm_ec*4.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-2);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-1);

                IntVect ip(ic); ip[0] += 1;
                int pR = (int) pmap_E_rhs[dit][dir](ip, 0);
                if(pR >= 0) {
                  vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane-2);
                }
              }
              else if(m_use_filtering && n==3) { // dJ_{ic}/dE_{ic-3}
                IntVect im(ic); im[0] -= 1;
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane-2);
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane-1);

                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-2);
              }
              else if(m_use_filtering && n==4) { // dJ_{ic}/dE_{ic-4}
                IntVect im(ic); im[0] -= 1;
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane-2);
              }
              else {
                vals[ix] = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-n);
              }
#else
              vals[ix] = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-n);
#endif
              ix++;
            }

            IntVect iR(ic); iR[0] += n;
            int pR = (int) pmap_E_rhs[dit][dir](iR, 0);
            if (pR >= 0) {
              icols[ix] = pR;
              vals[ix] = 0.0;
#ifdef FILTER_PC_TEST
              if(m_use_filtering && n==1) { // dJ_{ic}/dE_{ic+1}
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);
                vals[ix] += sig_norm_ec*4.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+1);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+2);

                IntVect im(ic); im[0] -= 1;
                int pL = (int) pmap_E_rhs[dit][dir](im, 0);
                if(pL >= 0) {
                  vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane+1);
                  vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane+2);
                }
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iR,diag_comp_inPlane-1);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](iR,diag_comp_inPlane);
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](iR,diag_comp_inPlane+1);
              }
              else if(m_use_filtering && n==2) { // dJ_{ic}/dE_{ic+2}
                IntVect im(ic); im[0] -= 1;
                int pL = (int) pmap_E_rhs[dit][dir](im, 0);
                if(pL >= 0) {
                  vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](im,diag_comp_inPlane+2);
                }

                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+1);
                vals[ix] += sig_norm_ec*4.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+2);

                IntVect ip(ic); ip[0] += 1;
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane+1);
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane+2);
              }
              else if(m_use_filtering && n==3) { // dJ_{ic}/dE_{ic+3}
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+2);

                IntVect ip(ic); ip[0] += 1;
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane+1);
                vals[ix] += sig_norm_ec*2.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane+2);
              }
              else if(m_use_filtering && n==4) { // dJ_{ic}/dE_{ic+4}
                IntVect ip(ic); ip[0] += 1;
                vals[ix] += sig_norm_ec*1.0/16.0*m_sigma_xx_pc[dit][dir](ip,diag_comp_inPlane+2);
              }
              else {
                vals[ix] = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+n);
              }
#else
              vals[ix] = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+n);
#endif
              ix++;
            }
          }

          if (m_pc_mass_matrix_include_ij) {

            /* sigma_xy contributions */
            for (int n = 0; n < m_pc_mass_matrix_width; n++) {

              if (n > m_sigma_xy_pc.nComp()/2-1) break;

              IntVect iL(ic); iL[0] -= n;
              int pL = (int) pmap_Ev_rhs[dit](iL,0);
              if (pL >= 0) {
                icols[ix] = pL;
                vals[ix] = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,m_sigma_xy_pc.nComp()/2-1-n);
                ix++;
              }

              IntVect iR(ic); iR[0] += (n+1);
              int pR = (int) pmap_Ev_rhs[dit](iR,0);
              if (pR >= 0) {
                icols[ix] = pR;
                vals[ix] = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,m_sigma_xy_pc.nComp()/2+n);
                ix++;
              }

            }

            /* sigma_xz contributions */
            for (int n = 0; n < m_pc_mass_matrix_width; n++) {

              if (n > m_sigma_xz_pc.nComp()/2-1) break;

              IntVect iL(ic); iL[0] -= n;
              int pL = (int) pmap_Ev_rhs[dit](iL,1);
              if (pL >= 0) {
                icols[ix] = pL;
                vals[ix] = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,m_sigma_xz_pc.nComp()/2-1-n);
                ix++;
              }

              IntVect iR(ic); iR[0] += (n+1);
              int pR = (int) pmap_Ev_rhs[dit](iR,1);
              if (pR >= 0) {
                icols[ix] = pR;
                vals[ix] = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,m_sigma_xz_pc.nComp()/2+n);
                ix++;
              }

            }
          }
        }
#elif CH_SPACEDIM==2
        if (    a_use_mass_matrices
             && (    (m_advanceE_comp[0] && (dir==0))
                  || (m_advanceE_comp[1] && (dir==1)) ) ) {

          /* sigma_xx/yy contributions */
          int idx(0);
          int pc_mass_matrix_width_x = 0;
          int pc_mass_matrix_width_y = 0;
          if (dir == 0) {
            pc_mass_matrix_width_x = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:2);
            pc_mass_matrix_width_y = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:3);
          } else {
            pc_mass_matrix_width_x = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:3);
            pc_mass_matrix_width_y = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:2);
          }
          for (int j = 0; j < (1+2*pc_mass_matrix_width_y); j++) {
            for (int i = 0; i < (1+2*pc_mass_matrix_width_x); i++) {
              if (!((i==pc_mass_matrix_width_x) && (j==pc_mass_matrix_width_y))) {
                IntVect iL(ic);
                iL[0] += (i-pc_mass_matrix_width_x);
                iL[1] += (j-pc_mass_matrix_width_y);
                int pL = (int) pmap_E_rhs[dit][dir](iL, 0);
                if (pL >= 0) {
                  icols[ix] = pL;
                  vals[ix] = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,idx);
                  ix++;
                }
              }
              idx++;
            }
          }

          if (m_pc_mass_matrix_include_ij) {

            if (dir == 0) {

              /* sigma_xy contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>3?3:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<2*pc_mass_matrix_width; j++) {
                  for (int i(0); i<2*pc_mass_matrix_width; i++) {
                    IntVect iL(ic);
                    iL[0] += (i+1-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_E_rhs[dit][dir+1](iL,0);
                    if (pL >= 0) {
                      icols[ix] = pL;
                      vals[ix] = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,idx);
                    }
                    idx++;
                  }
                }
              }
              /* sigma_xz contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>2?2:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<(1+2*pc_mass_matrix_width); j++) {
                  for (int i(0); i<2*pc_mass_matrix_width; i++) {
                    IntVect iL(ic);
                    iL[0] += (i+1-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_Ev_rhs[dit](iL,0);
                    if (pL >= 0) {
                      icols[ix] = pL;
                      vals[ix] = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,idx);
                    }
                    idx++;
                  }
                }
              }

            } else {

              /* sigma_yx contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>3?3:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<2*pc_mass_matrix_width; j++) {
                  for (int i(0); i<2*pc_mass_matrix_width; i++) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j+1-pc_mass_matrix_width);
                    int pL = (int) pmap_E_rhs[dit][dir-1](iL,0);
                    if (pL >= 0) {
                      icols[ix] = pL;
                      vals[ix] = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,idx);
                    }
                    idx++;
                  }
                }
              }
              /* sigma_yz contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>2?2:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<2*pc_mass_matrix_width; j++) {
                  for (int i(0); i<(1+2*pc_mass_matrix_width); i++) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j+1-pc_mass_matrix_width);
                    int pL = (int) pmap_Ev_rhs[dit](iL,0);
                    if (pL >= 0) {
                      icols[ix] = pL;
                      vals[ix] = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,idx);
                    }
                    idx++;
                  }
                }
              }

            }
          }
        }
#endif

#if CH_SPACEDIM==2
        if (a_include_EM) {
          int tdir = (dir + 1) % SpaceDim;
          int sign = 1 - 2*dir;
          if(m_mesh.anticyclic()) sign *= -1;

	      // JRA, here is where axisymmetric modifications are needed in 2D RZ
          // dEz_{i,j+1/2}/dt + Jz_{i,j+1/2} = (r_{i+1/2}*Bth_{i+1/2,j+1/2} - r_{i-1/2}*Bth_{i-1/2,j+1/2})/dr/r_{i}
          // for r_{i=0}=0, d(rBth)/dr/r = 4*Bth_{i=1/2,j+1/2}/dr

          IntVect iL(ic); iL[tdir]--;
          Real aL = (-1) * sign * (cnormDt/dX[tdir]);
          aL *= (m_advanceE_comp[dir] * m_PC_mask_Bv[dit](iL,2*tdir));
          int pL = (int) pmap_Bv_rhs[dit](iL, 0);

          if (pL >= 0) {
            icols[ix] = pL;
            if(m_mesh.axisymmetric() && dir==1) {
              Real local_Xec = Xec[dit][dir](ic,0);
              if(local_Xec==0.0) vals[ix] = 0.0;
              else vals[ix] = -aL*Xcc[dit](iL,0)/local_Xec;
            }
            else vals[ix] = -aL;
            ix++;
          }

          IntVect iR(ic);
          Real aR = sign * (cnormDt/dX[tdir]);
          aR *=  (m_advanceE_comp[dir] * m_PC_mask_Bv[dit](iR,2*tdir+1));
          int pR = (int) pmap_Bv_rhs[dit](iR, 0);

          if (pR >= 0) {
            icols[ix] = pR;
            if(m_mesh.axisymmetric() && dir==1) {
              Real local_Xec = Xec[dit][dir](ic,0);
              if(local_Xec==0.0) vals[ix] = -aR*4.0;
              else vals[ix] = -aR*Xcc[dit](iR,0)/local_Xec;
            }
            else vals[ix] = -aR;
            ix++;
          }
        }
#endif

        CH_assert(ix <= m_pcmat_nnz);
        CH_assert(ix <= a_P.getNBands());

        a_P.setRowValues( pc, ix, icols.data(), vals.data() );
      }
    }

#if CH_SPACEDIM<3
    /* virtual electric field */
    {
      const NodeFArrayBox& pmap( pmap_Ev_lhs[dit] );

      auto box( surroundingNodes( grow(pmap.box(), -ghosts) ) );
      for (int dir=0; dir<SpaceDim; ++dir) {
         int idir_bdry = phys_domain.domainBox().bigEnd(dir);
         if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {

        auto ic = bit();
        const Real sig_norm_nc = sig_normC/Jnc[dit](ic,0);

        for (int n(0); n < pmap.nComp(); n++) {

          int pc = (int) pmap(ic, n);

          int ncols = m_pcmat_nnz;
          std::vector<int> icols(ncols, -1);
          std::vector<Real> vals(ncols, 0.0);
          int ix = 0;

          icols[ix] = pc;
          vals[ix] = 1.0;
#if CH_SPACEDIM==1
          if(a_use_mass_matrices && m_advanceE_comp[1+n]) {
              vals[ix] += sig_norm_nc * ( n==0 ? m_sigma_yy_pc[dit](ic,diag_comp_virtual)
                                               : m_sigma_zz_pc[dit](ic,diag_comp_virtual) );
          }
#elif CH_SPACEDIM==2
          if(a_use_mass_matrices && m_advanceE_comp[2]) {
            vals[ix] += sig_norm_nc*m_sigma_zz_pc[dit](ic,diag_comp_virtual);
          }
#endif
          ix++;

          if(a_use_mass_matrices && m_advanceE_comp[SpaceDim+n]) {
#if CH_SPACEDIM==1
            /* sigma_yy/zz contributions */
            for (int s = 1; s < m_pc_mass_matrix_width+1; s++) {

              if (s > (m_sigma_yy_pc.nComp()-1)/2) break;
              if (s > (m_sigma_zz_pc.nComp()-1)/2) break;

              IntVect iL(ic); iL[0] -= s;
              int pL = (int) pmap_Ev_rhs[dit](iL, n);
              if (pL >= 0) {
                icols[ix] = pL;
                vals[ix] = sig_norm_nc * ( n==0 ? m_sigma_yy_pc[dit](ic,diag_comp_virtual-s)
                                                : m_sigma_zz_pc[dit](ic,diag_comp_virtual-s) );
                ix++;
              }

              IntVect iR(ic); iR[0] += s;
              int pR = (int) pmap_Ev_rhs[dit](iR, n);
              if (pR >= 0) {
                icols[ix] = pR;
                vals[ix] = sig_norm_nc * ( n==0 ? m_sigma_yy_pc[dit](ic,diag_comp_virtual+s)
                                                : m_sigma_zz_pc[dit](ic,diag_comp_virtual+s) );
                ix++;
              }
            }
            if (m_pc_mass_matrix_include_ij) {

              // sigma_yx/zx contributions
              for (int s = 0; s < m_pc_mass_matrix_width; s++) {

                if (s > m_sigma_yx_pc.nComp()/2-1) break;

                IntVect iL(ic); iL[0] -= (s+1);
                int pL = (int) pmap_E_rhs[dit][0](iL,0);
                if (pL >= 0) {
                  icols[ix] = pL;
                  vals[ix] = sig_norm_nc * ( n==0 ? m_sigma_yx_pc[dit](ic,m_sigma_yx_pc.nComp()/2-1-s)
                                                  : m_sigma_zx_pc[dit](ic,m_sigma_zx_pc.nComp()/2-1-s) );
                  ix++;
                }

                IntVect iR(ic); iR[0] += s;
                int pR = (int) pmap_E_rhs[dit][0](iR,0);
                if (pR >= 0) {
                  icols[ix] = pR;
                  vals[ix] = sig_norm_nc * ( n==0 ? m_sigma_yx_pc[dit](ic,m_sigma_yx_pc.nComp()/2+s)
                                                  : m_sigma_zx_pc[dit](ic,m_sigma_zx_pc.nComp()/2+s) );
                  ix++;
                }

              }

              // sigma_yz/zy contributions
              {
                IntVect ii(ic);
                int pI = (int) pmap_Ev_rhs[dit](ii, !n);
                if (pI >= 0) {
                  icols[ix] = pI;
                  vals[ix] = sig_norm_nc * ( n==0 ? m_sigma_yz_pc[dit](ic,(m_sigma_yz_pc.nComp()-1)/2)
                                                  : m_sigma_zy_pc[dit](ic,(m_sigma_zy_pc.nComp()-1)/2) );
                  ix++;
                }
              }
              for (int s = 1; s < m_pc_mass_matrix_width+1; s++) {

                if (s > (m_sigma_yz_pc.nComp()-1)/2) break;

                IntVect iL(ic); iL[0] -= s;
                int pL = (int) pmap_Ev_rhs[dit](iL, !n);
                if (pL >= 0) {
                  icols[ix] = pL;
                  vals[ix] = sig_norm_nc * ( n==0 ? m_sigma_yz_pc[dit](ic,(m_sigma_yz_pc.nComp()-1)/2-s)
                                                  : m_sigma_zy_pc[dit](ic,(m_sigma_zy_pc.nComp()-1)/2-s) );
                  ix++;
                }

                IntVect iR(ic); iR[0] += s;
                int pR = (int) pmap_Ev_rhs[dit](iR, !n);
                if (pR >= 0) {
                  icols[ix] = pR;
                  vals[ix] = sig_norm_nc * ( n==0 ? m_sigma_yz_pc[dit](ic,(m_sigma_yz_pc.nComp()-1)/2+s)
                                                  : m_sigma_zy_pc[dit](ic,(m_sigma_zy_pc.nComp()-1)/2+s) );
                  ix++;
                }
              }

            }
#elif CH_SPACEDIM==2

            /* sigma_zz contributions */
            {
              int idx(0);
              int pc_mass_matrix_width(m_pc_mass_matrix_width<2?m_pc_mass_matrix_width:1);
              for (int j = 0; j < (1+2*pc_mass_matrix_width); j++) {
                for (int i = 0; i < (1+2*pc_mass_matrix_width); i++) {
                  if (!((i==pc_mass_matrix_width) && (j==pc_mass_matrix_width))) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_Ev_rhs[dit](iL,n);
                    if (pL >= 0) {
                      icols[ix] = pL;
                      vals[ix] = sig_norm_nc *  m_sigma_zz_pc[dit](ic,idx);
                      ix++;
                    }
                  }
                  idx++;
                }
              }
            }

            if (m_pc_mass_matrix_include_ij) {
              /* sigma_zx/zy contribution */
              for (int dir = 0; dir < SpaceDim; dir++) {
                int idx(0);
                int pc_mass_matrix_width(m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:2);
                for (int j(0); j<(dir+2*pc_mass_matrix_width); j++) {
                  for (int i(0); i<((1-dir)+2*pc_mass_matrix_width); i++) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_E_rhs[dit][dir](iL,0);
                    if (pL >= 0) {
                      icols[ix] = pL;
                      vals[ix] = sig_norm_nc * (dir == 0 ? m_sigma_zx_pc[dit](ic,idx)
                                                         : m_sigma_zy_pc[dit](ic,idx) );
                    }
                    idx++;
                  }
                }
              }
            }
#endif
          }

          if (a_include_EM) {
#if CH_SPACEDIM==1

            int dirj = 0, dirk = 1;
            if (n == 0) {

              int sign = -1;

              // axisymmetric modifications are needed in 1D sph geom
              // sph: dEth_{i}/dt + Jth_{i} = -(r_{i+1/2}*Bphi_{i+1/2} - r_{i-1/2}*Bphi_{i-1/2})/dr/r_{i}
              // for r_{i=0}=0, d(rBphi)/dr/r = 4*Bphi_{i=1/2}/dr

              Real aL = (-1) * sign * (cnormDt/dX[0]) * m_advanceE_comp[n+1];
              IntVect iL(ic); iL[0]--;
              int pL = (int) pmap_Bv_rhs[dit](iL, dirk);

              if (pL >= 0) {
                icols[ix] = pL;
                if(m_mesh.axisymmetric() && geom_type=="sph_R") {
                  Real local_Xnc = Xnc[dit](ic,0);
                  if(local_Xnc==0.0) { vals[ix] = 0.0; }
                  else { vals[ix] = -aL*Xcc[dit](iL,0)/local_Xnc; }
                }
                else { vals[ix] = -aL; }
                ix++;
              }

              Real aR =  sign * (cnormDt/dX[0]) * m_advanceE_comp[n+1];
              IntVect iR(ic);
              int pR = (int) pmap_Bv_rhs[dit](iR, dirk);

              if (pR >= 0) {
                icols[ix] = pR;
                if(m_mesh.axisymmetric() && geom_type=="sph_R") {
                  Real local_Xnc = Xnc[dit](ic,0);
                  if(local_Xnc==0.0) { vals[ix] = -aR*4.0; }
                  else { vals[ix] = -aR*Xcc[dit](iR,0)/local_Xnc; }
                }
                else { vals[ix] = -aR; }
                ix++;
              }

            } else if (n == 1) {

              int sign = 1;

	          // axisymmetric modifications are needed for 1D cyl and 1D sph geoms
              // cyl: dEz_{i}/dt + Jz_{i} = (r_{i+1/2}*Bth_{i+1/2} - r_{i-1/2}*Bth_{i-1/2})/dr/r_{i}
              // sph: dEphi_{i}/dt + Jphi_{i} = (r_{i+1/2}*Bth_{i+1/2} - r_{i-1/2}*Bth_{i-1/2})/dr/r_{i}
              // for r_{i=0}=0, d(rBth)/dr/r = 4*Bth_{i=1/2}/dr

              Real aL = (-1) * sign * (cnormDt/dX[0]) * m_advanceE_comp[n+1];
              IntVect iL(ic); iL[0]--;
              int pL = (int) pmap_Bv_rhs[dit](iL, dirj);

              if (pL >= 0) {
                icols[ix] = pL;
                if(m_mesh.axisymmetric()) {
                  Real local_Xnc = Xnc[dit](ic,0);
                  if(local_Xnc==0.0) { vals[ix] = 0.0; }
                  else { vals[ix] = -aL*Xcc[dit](iL,0)/local_Xnc; }
                }
                else { vals[ix] = -aL; }
                ix++;
              }

              Real aR =  sign * (cnormDt/dX[0]) * m_advanceE_comp[n+1];
              IntVect iR(ic);
              int pR = (int) pmap_Bv_rhs[dit](iR, dirj);

              if (pR >= 0) {
                icols[ix] = pR;
                if(m_mesh.axisymmetric()) {
                  Real local_Xnc = Xnc[dit](ic,0);
                  if(local_Xnc==0.0) { vals[ix] = -aR*4.0; }
                  else { vals[ix] = -aR*Xcc[dit](iR,0)/local_Xnc; }
                }
                else { vals[ix] = -aR; }
                ix++;
              }

            }
#elif CH_SPACEDIM==2
            int dirj = 0, dirk = 1;
            {
              int sign = 1;
              if(m_mesh.anticyclic()) sign *= -1;

              IntVect iL(ic); iL[dirj]--;
              Real aL = (-1) * sign * (cnormDt/dX[dirj]);
              aL *= (m_advanceE_comp[2] * m_PC_mask_B[dit][dirk](iL,0));
              iL[dirj] += m_PC_mask_B[dit][dirk](iL,m_PC_mask_B.nComp()-2);
              int pL = (int) pmap_B_rhs[dit][dirk](iL, 0);

              if (pL >= 0) {
                icols[ix] = pL;
                vals[ix] = -aL;
                ix++;
              }

              IntVect iR(ic);
              Real aR = sign * (cnormDt/dX[dirj]);
              aR *= (m_advanceE_comp[2] * m_PC_mask_B[dit][dirk](iR,1));
              iR[dirj] += m_PC_mask_B[dit][dirk](iR,m_PC_mask_B.nComp()-1);
              int pR = (int) pmap_B_rhs[dit][dirk](iR, 0);

              if (pR >= 0) {
                icols[ix] = pR;
                vals[ix] = -aR;
                ix++;
              }
            }
            {
              int sign = -1;
              if(m_mesh.anticyclic()) sign *= -1;

              IntVect iL(ic); iL[dirk]--;
              Real aL = (-1) * sign * (cnormDt/dX[dirk]);
              aL *= (m_PC_mask_B[dit][dirj](iL,0) * m_advanceE_comp[2]);
              int pL = (int) pmap_B_rhs[dit][dirj](iL, 0);

              if (pL >= 0) {
                icols[ix] = pL;
                vals[ix] = -aL;
                ix++;
              }

              IntVect iR(ic);
              Real aR = sign * (cnormDt/dX[dirk]);
              aR *= (m_PC_mask_B[dit][dirj](iR,1) * m_advanceE_comp[2]);
              int pR = (int) pmap_B_rhs[dit][dirj](iR, 0);

              if (pR >= 0) {
                icols[ix] = pR;
                vals[ix] = -aR;
                ix++;
              }
            }
#endif
          }

          CH_assert(ix <= m_pcmat_nnz);
          CH_assert(ix <= a_P.getNBands());

          a_P.setRowValues( pc, ix, icols.data(), vals.data() );
        }
      }
    }
#endif

  }

  return;
}

void EMFields::assemblePrecondMatrixCurl2(  BandedMatrix&  a_P,
                                      const bool           a_include_EM,
                                      const bool           a_use_mass_matrices,
                                      const Real           a_dt )
{
  CH_TIME("EMFields::assemblePrecondMatrixE()");

  auto grids(m_mesh.getDBL());

  // copy over the mass matrix elements used in the PC and do addOp exchange()
  const Real cnormDt = (a_dt*m_cvacNorm) * (a_dt*m_cvacNorm); //TODO ask JRA
  const Real sig_normC = (a_dt*m_cvacNorm)*m_Jnorm_factor;
  int diag_comp_inPlane=0, diag_comp_virtual=0;
  if(a_use_mass_matrices) {
    diag_comp_inPlane = (m_sigma_xx_pc.nComp()-1)/2;
#if CH_SPACEDIM<3
    diag_comp_virtual = (m_sigma_zz_pc.nComp()-1)/2;
#endif
  }

  const LevelData<EdgeDataBox>& pmap_E_lhs = m_gdofs.dofDataLHSE();
  const LevelData<NodeFArrayBox>& pmap_Ev_lhs = m_gdofs.dofDataLHSEv();
  const LevelData<EdgeDataBox>& pmap_E_rhs = m_gdofs.dofDataRHSE();
  const LevelData<NodeFArrayBox>& pmap_Ev_rhs = m_gdofs.dofDataRHSEv();

#if CH_SPACEDIM<3
  const LevelData<NodeFArrayBox>& Jnc = m_mesh.getCorrectedJnc();
#endif
  const LevelData<EdgeDataBox>& Jec = m_mesh.getCorrectedJec();

  auto phys_domain( grids.physDomain() );
  const int ghosts(m_mesh.ghosts());
  const RealVect& dX(m_mesh.getdX());


#if CH_SPACEDIM==1
  const string& geom_type = m_mesh.geomType();
#endif

  for (auto dit( grids.dataIterator() ); dit.ok(); ++dit) {

    /* electric field */
    for (int dir = 0; dir < SpaceDim; dir++) {

      const FArrayBox& pmap( pmap_E_lhs[dit][dir] );

      auto box( grow(pmap.box(), -ghosts) );
      for (int adir=0; adir<SpaceDim; ++adir) {
         if (adir != dir) {
            int idir_bdry = phys_domain.domainBox().bigEnd(adir);
            if (box.bigEnd(adir) < idir_bdry) box.growHi(adir, -1);
         }
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {

        auto ic = bit();
        CH_assert( pmap.nComp() == 1);
        int pc = (int) pmap(ic, 0);

	    const Real sig_norm_ec = sig_normC/Jec[dit][dir](ic,0);
        std::map<int,Real> this_row;

        {
          Real val = 1.0;
          if (a_use_mass_matrices) { // dJ_{ic}/dE_{ic}
             val += sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane);
          }
          insertOrAdd<int,Real>(this_row, pc, val);
        }

        /* curl-curl terms */
        if (a_include_EM) {
#if CH_SPACEDIM==2
          {
            IntVect ij(ic);
            Real aj = 0.0;
            if (dir == 0) {
              aj = 2*cnormDt/(dX[1]*dX[1])*m_advanceE_comp[dir]*m_PC_mask_E[dit][dir](ij,0);
            } else if (dir == 1) {
              aj = 2*cnormDt/(dX[0]*dX[0])*m_advanceE_comp[dir]*m_PC_mask_E[dit][dir](ij,0);
            }
            aj *= m_PC_mask_E[dit][dir](ij,0);
            int pj = (int) pmap_E_rhs[dit][dir](ij, 0);
            if (pj >= 0) {
              insertOrAdd<int,Real>(this_row, pj, aj);
            }
          }
          if (dir < 2) {
            int tdir = (dir == 0 ? 1 : 0);
            {
              IntVect ij(ic); ij[tdir]--;
              Real aj = -cnormDt/(dX[tdir]*dX[tdir]) * (m_advanceE_comp[dir]);
              aj *= m_PC_mask_E[dit][dir](ij,2*tdir);
              aj *= m_PC_mask_E[dit][dir](ic,1);
              int pj = (int) pmap_E_rhs[dit][dir](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
            {
              IntVect ij(ic); ij[tdir]++;
              Real aj = -cnormDt/(dX[tdir]*dX[tdir]) * (m_advanceE_comp[dir]);
              aj *= m_PC_mask_E[dit][dir](ij,2*tdir+1);
              aj *= m_PC_mask_E[dit][dir](ic,1);
              int pj = (int) pmap_E_rhs[dit][dir](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
            {
              IntVect ij(ic); ij[tdir]--;
              Real aj = cnormDt/(dX[dir]*dX[tdir]) * (m_advanceE_comp[dir]);
              aj *= m_PC_mask_E[dit][tdir](ij,2*tdir);
              aj *= m_PC_mask_E[dit][dir](ic,1);
              int pj = (int) pmap_E_rhs[dit][tdir](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
            {
              IntVect ij(ic);
              Real aj = -cnormDt/(dX[dir]*dX[tdir]) * (m_advanceE_comp[dir]);
              aj *= m_PC_mask_E[dit][tdir](ij,2*tdir);
              aj *= m_PC_mask_E[dit][dir](ic,1);
              int pj = (int) pmap_E_rhs[dit][tdir](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
            {
              IntVect ij(ic); ij[dir]++; ij[tdir]--;
              Real aj = -cnormDt/(dX[dir]*dX[tdir]) * (m_advanceE_comp[dir]);
              aj *= m_PC_mask_E[dit][tdir](ij,2*tdir+1);
              aj *= m_PC_mask_E[dit][dir](ic,1);
              int pj = (int) pmap_E_rhs[dit][tdir](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
            {
              IntVect ij(ic); ij[dir]++;;
              Real aj = cnormDt/(dX[dir]*dX[tdir]) * (m_advanceE_comp[dir]);
              aj *= m_PC_mask_E[dit][tdir](ij,2*tdir+1);
              aj *= m_PC_mask_E[dit][dir](ic,1);
              int pj = (int) pmap_E_rhs[dit][tdir](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
          }
#endif
        }

#if CH_SPACEDIM==1
        if (a_use_mass_matrices && m_advanceE_comp[0]) {

          /* sigma_xx contributions */
          for (int n = 1; n < m_pc_mass_matrix_width+1; n++) {

            if (n > (m_sigma_xx_pc.nComp()-1)/2) break;

            IntVect iL(ic); iL[0] -= n;
            int pL = (int) pmap_E_rhs[dit][dir](iL, 0);
            if (pL >= 0) {
              auto val = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane-n);
              insertOrAdd<int,Real>(this_row, pL, val);
            }

            IntVect iR(ic); iR[0] += n;
            int pR = (int) pmap_E_rhs[dit][dir](iR, 0);
            if (pR >= 0) {
              auto val = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,diag_comp_inPlane+n);
              insertOrAdd<int,Real>(this_row, pR, val);
            }
          }

          if (m_pc_mass_matrix_include_ij) {

            /* sigma_xy contributions */
            for (int n = 0; n < m_pc_mass_matrix_width; n++) {

              if (n > m_sigma_xy_pc.nComp()/2-1) break;

              IntVect iL(ic); iL[0] -= n;
              int pL = (int) pmap_Ev_rhs[dit](iL,0);
              if (pL >= 0) {
                auto val = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,m_sigma_xy_pc.nComp()/2-1-n);
                insertOrAdd<int,Real>(this_row, pL, val);
              }

              IntVect iR(ic); iR[0] += (n+1);
              int pR = (int) pmap_Ev_rhs[dit](iR,0);
              if (pR >= 0) {
                auto val = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,m_sigma_xy_pc.nComp()/2+n);
                insertOrAdd<int,Real>(this_row, pR, val);
              }

            }

            /* sigma_xz contributions */
            for (int n = 0; n < m_pc_mass_matrix_width; n++) {

              if (n > m_sigma_xz_pc.nComp()/2-1) break;

              IntVect iL(ic); iL[0] -= n;
              int pL = (int) pmap_Ev_rhs[dit](iL,1);
              if (pL >= 0) {
                auto val = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,m_sigma_xz_pc.nComp()/2-1-n);
                insertOrAdd<int,Real>(this_row, pL, val);
              }

              IntVect iR(ic); iR[0] += (n+1);
              int pR = (int) pmap_Ev_rhs[dit](iR,1);
              if (pR >= 0) {
                auto val = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,m_sigma_xz_pc.nComp()/2+n);
                insertOrAdd<int,Real>(this_row, pR, val);
              }

            }
          }
        }
#elif CH_SPACEDIM==2
        if (    a_use_mass_matrices
             && (    (m_advanceE_comp[0] && (dir==0))
                  || (m_advanceE_comp[1] && (dir==1)) ) ) {

          /* sigma_xx/yy contributions */
          int idx(0);
          int pc_mass_matrix_width_x = 0;
          int pc_mass_matrix_width_y = 0;
          if (dir == 0) {
            pc_mass_matrix_width_x = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:2);
            pc_mass_matrix_width_y = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:3);
          } else {
            pc_mass_matrix_width_x = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:3);
            pc_mass_matrix_width_y = (m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:2);
          }
          for (int j = 0; j < (1+2*pc_mass_matrix_width_y); j++) {
            for (int i = 0; i < (1+2*pc_mass_matrix_width_x); i++) {
              if (!((i==pc_mass_matrix_width_x) && (j==pc_mass_matrix_width_y))) {
                IntVect iL(ic);
                iL[0] += (i-pc_mass_matrix_width_x);
                iL[1] += (j-pc_mass_matrix_width_y);
                int pL = (int) pmap_E_rhs[dit][dir](iL, 0);
                if (pL >= 0) {
                  auto val = sig_norm_ec*m_sigma_xx_pc[dit][dir](ic,idx);
                  insertOrAdd<int,Real>(this_row, pL, val);
                }
              }
              idx++;
            }
          }

          if (m_pc_mass_matrix_include_ij) {

            if (dir == 0) {

              /* sigma_xy contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>3?3:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<2*pc_mass_matrix_width; j++) {
                  for (int i(0); i<2*pc_mass_matrix_width; i++) {
                    IntVect iL(ic);
                    iL[0] += (i+1-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_E_rhs[dit][dir+1](iL,0);
                    if (pL >= 0) {
                      auto val = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,idx);
                      insertOrAdd<int,Real>(this_row, pL, val);
                    }
                    idx++;
                  }
                }
              }
              /* sigma_xz contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>2?2:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<(1+2*pc_mass_matrix_width); j++) {
                  for (int i(0); i<2*pc_mass_matrix_width; i++) {
                    IntVect iL(ic);
                    iL[0] += (i+1-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_Ev_rhs[dit](iL,0);
                    if (pL >= 0) {
                      auto val = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,idx);
                      insertOrAdd<int,Real>(this_row, pL, val);
                    }
                    idx++;
                  }
                }
              }

            } else {

              /* sigma_yx contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>3?3:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<2*pc_mass_matrix_width; j++) {
                  for (int i(0); i<2*pc_mass_matrix_width; i++) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j+1-pc_mass_matrix_width);
                    int pL = (int) pmap_E_rhs[dit][dir-1](iL,0);
                    if (pL >= 0) {
                      auto val = sig_norm_ec*m_sigma_xy_pc[dit][dir](ic,idx);
                      insertOrAdd<int,Real>(this_row, pL, val);
                    }
                    idx++;
                  }
                }
              }
              /* sigma_yz contributions */
              {
                int pc_mass_matrix_width(m_pc_mass_matrix_width>2?2:m_pc_mass_matrix_width);
                int idx(0);
                for (int j(0); j<2*pc_mass_matrix_width; j++) {
                  for (int i(0); i<(1+2*pc_mass_matrix_width); i++) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j+1-pc_mass_matrix_width);
                    int pL = (int) pmap_Ev_rhs[dit](iL,0);
                    if (pL >= 0) {
                      auto val = sig_norm_ec*m_sigma_xz_pc[dit][dir](ic,idx);
                      insertOrAdd<int,Real>(this_row, pL, val);
                    }
                    idx++;
                  }
                }
              }

            }
          }
        }
#endif

        CH_assert(this_row.size() <= m_pcmat_nnz);
        CH_assert(this_row.size() <= a_P.getNBands());

        std::vector<int> icols(0);
        std::vector<Real> vals(0);
        for (auto it = this_row.cbegin(); it != this_row.cend(); ++it) {
            icols.push_back( it->first );
            vals.push_back( it->second );
        }

        a_P.setRowValues( pc, this_row.size(), icols.data(), vals.data() );
      }
    }

#if CH_SPACEDIM<3
    /* virtual electric field */
    {
      const NodeFArrayBox& pmap( pmap_Ev_lhs[dit] );

      auto box( surroundingNodes( grow(pmap.box(), -ghosts) ) );
      for (int dir=0; dir<SpaceDim; ++dir) {
         int idir_bdry = phys_domain.domainBox().bigEnd(dir);
         if (box.bigEnd(dir) < idir_bdry) box.growHi(dir, -1);
      }

      for (BoxIterator bit(box); bit.ok(); ++bit) {

        auto ic = bit();
        const Real sig_norm_nc = sig_normC/Jnc[dit](ic,0);

        for (int n(0); n < pmap.nComp(); n++) {

          std::map<int,Real> this_row;
          int pc = (int) pmap(ic, n);

          {
            Real val = 1.0;
#if CH_SPACEDIM==1
            if(a_use_mass_matrices && m_advanceE_comp[1+n]) {
                val += sig_norm_nc * ( n==0 ? m_sigma_yy_pc[dit](ic,diag_comp_virtual)
                                            : m_sigma_zz_pc[dit](ic,diag_comp_virtual) );
            }
#elif CH_SPACEDIM==2
            if(a_use_mass_matrices && m_advanceE_comp[2]) {
              val += sig_norm_nc*m_sigma_zz_pc[dit](ic,diag_comp_virtual);
            }
#endif
            insertOrAdd<int,Real>(this_row, pc, val);
          }

        /* curl-curl terms */
        if (a_include_EM) {
#if CH_SPACEDIM==1
          {
            IntVect ij(ic);
            Real aj = 2*cnormDt/(dX[0]*dX[0])*m_advanceE_comp[1+n];
            int pj = (int) pmap_Ev_rhs[dit](ij, n);
            if (pj >= 0) {
              insertOrAdd<int,Real>(this_row, pj, aj);
            }
          }
          int dir = n+1;
          {
            IntVect ij(ic); ij[0]--;
            Real aj = -cnormDt/(dX[0]*dX[0]) * (m_advanceE_comp[dir]);
            int pj = (int) pmap_Ev_rhs[dit](ij, n);
            if (pj >= 0) {
              insertOrAdd<int,Real>(this_row, pj, aj);
            }
          }
          {
            IntVect ij(ic); ij[0]++;
            Real aj = -cnormDt/(dX[0]*dX[0]) * (m_advanceE_comp[dir]);
            int pj = (int) pmap_Ev_rhs[dit](ij, n);
            if (pj >= 0) {
              insertOrAdd<int,Real>(this_row, pj, aj);
            }
          }
#elif CH_SPACEDIM==2
          {
            IntVect ij(ic);
            Real aj = 2*cnormDt*(1/(dX[0]*dX[0])+1/(dX[1]*dX[1]))*m_advanceE_comp[2];
            aj *= m_PC_mask_Ev[dit](ij,0);
            int pj = (int) pmap_Ev_rhs[dit](ij, 0);
            if (pj >= 0) {
              insertOrAdd<int,Real>(this_row, pj, aj);
            }
          }
          for (int dir=0; dir<SpaceDim; dir++) {
            {
              IntVect ij(ic); ij[dir]--;
              Real aj = -cnormDt/(dX[dir]*dX[dir]) * (m_advanceE_comp[2]);
              aj *= m_PC_mask_Ev[dit](ij,2+2*dir) * m_PC_mask_Ev[dit](ic,1);
              int pj = (int) pmap_Ev_rhs[dit](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
            {
              IntVect ij(ic); ij[dir]++;
              Real aj = -cnormDt/(dX[dir]*dX[dir]) * (m_advanceE_comp[2]);
              aj *= m_PC_mask_Ev[dit](ij,2+2*dir+1) * m_PC_mask_Ev[dit](ic,1);
              int pj = (int) pmap_Ev_rhs[dit](ij, 0);
              if (pj >= 0) {
                insertOrAdd<int,Real>(this_row, pj, aj);
              }
            }
          }
#endif
        }

          if(a_use_mass_matrices && m_advanceE_comp[SpaceDim+n]) {
#if CH_SPACEDIM==1
            /* sigma_yy/zz contributions */
            for (int s = 1; s < m_pc_mass_matrix_width+1; s++) {

              if (s > (m_sigma_yy_pc.nComp()-1)/2) break;
              if (s > (m_sigma_zz_pc.nComp()-1)/2) break;

              IntVect iL(ic); iL[0] -= s;
              int pL = (int) pmap_Ev_rhs[dit](iL, n);
              if (pL >= 0) {
                auto val = sig_norm_nc * ( n==0 ? m_sigma_yy_pc[dit](ic,diag_comp_virtual-s)
                                                : m_sigma_zz_pc[dit](ic,diag_comp_virtual-s) );
                insertOrAdd<int,Real>(this_row, pL, val);
              }

              IntVect iR(ic); iR[0] += s;
              int pR = (int) pmap_Ev_rhs[dit](iR, n);
              if (pR >= 0) {
                auto val = sig_norm_nc * ( n==0 ? m_sigma_yy_pc[dit](ic,diag_comp_virtual+s)
                                                : m_sigma_zz_pc[dit](ic,diag_comp_virtual+s) );
                insertOrAdd<int,Real>(this_row, pR, val);
              }
            }
            if (m_pc_mass_matrix_include_ij) {

              // sigma_yx/zx contributions
              for (int s = 0; s < m_pc_mass_matrix_width; s++) {

                if (s > m_sigma_yx_pc.nComp()/2-1) break;

                IntVect iL(ic); iL[0] -= (s+1);
                int pL = (int) pmap_E_rhs[dit][0](iL,0);
                if (pL >= 0) {
                  auto val = sig_norm_nc * ( n==0 ? m_sigma_yx_pc[dit](ic,m_sigma_yx_pc.nComp()/2-1-s)
                                                  : m_sigma_zx_pc[dit](ic,m_sigma_zx_pc.nComp()/2-1-s) );
                  insertOrAdd<int,Real>(this_row, pL, val);
                }

                IntVect iR(ic); iR[0] += s;
                int pR = (int) pmap_E_rhs[dit][0](iR,0);
                if (pR >= 0) {
                  auto val = sig_norm_nc * ( n==0 ? m_sigma_yx_pc[dit](ic,m_sigma_yx_pc.nComp()/2+s)
                                                  : m_sigma_zx_pc[dit](ic,m_sigma_zx_pc.nComp()/2+s) );
                  insertOrAdd<int,Real>(this_row, pR, val);
                }

              }

              // sigma_yz/zy contributions
              {
                IntVect ii(ic);
                int pI = (int) pmap_Ev_rhs[dit](ii, !n);
                if (pI >= 0) {
                  auto val = sig_norm_nc * ( n==0 ? m_sigma_yz_pc[dit](ic,(m_sigma_yz_pc.nComp()-1)/2)
                                                  : m_sigma_zy_pc[dit](ic,(m_sigma_zy_pc.nComp()-1)/2) );
                  insertOrAdd<int,Real>(this_row, pI, val);
                }
              }
              for (int s = 1; s < m_pc_mass_matrix_width+1; s++) {

                if (s > (m_sigma_yz_pc.nComp()-1)/2) break;

                IntVect iL(ic); iL[0] -= s;
                int pL = (int) pmap_Ev_rhs[dit](iL, !n);
                if (pL >= 0) {
                  auto val = sig_norm_nc * ( n==0 ? m_sigma_yz_pc[dit](ic,(m_sigma_yz_pc.nComp()-1)/2-s)
                                                  : m_sigma_zy_pc[dit](ic,(m_sigma_zy_pc.nComp()-1)/2-s) );
                  insertOrAdd<int,Real>(this_row, pL, val);
                }

                IntVect iR(ic); iR[0] += s;
                int pR = (int) pmap_Ev_rhs[dit](iR, !n);
                if (pR >= 0) {
                  auto val = sig_norm_nc * ( n==0 ? m_sigma_yz_pc[dit](ic,(m_sigma_yz_pc.nComp()-1)/2+s)
                                                  : m_sigma_zy_pc[dit](ic,(m_sigma_zy_pc.nComp()-1)/2+s) );
                  insertOrAdd<int,Real>(this_row, pR, val);
                }
              }

            }
#elif CH_SPACEDIM==2

            /* sigma_zz contributions */
            {
              int idx(0);
              int pc_mass_matrix_width(m_pc_mass_matrix_width<2?m_pc_mass_matrix_width:1);
              for (int j = 0; j < (1+2*pc_mass_matrix_width); j++) {
                for (int i = 0; i < (1+2*pc_mass_matrix_width); i++) {
                  if (!((i==pc_mass_matrix_width) && (j==pc_mass_matrix_width))) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_Ev_rhs[dit](iL,n);
                    if (pL >= 0) {
                      auto val = sig_norm_nc *  m_sigma_zz_pc[dit](ic,idx);
                      insertOrAdd<int,Real>(this_row, pL, val);
                    }
                  }
                  idx++;
                }
              }
            }

            if (m_pc_mass_matrix_include_ij) {
              /* sigma_zx/zy contribution */
              for (int dir = 0; dir < SpaceDim; dir++) {
                int idx(0);
                int pc_mass_matrix_width(m_pc_mass_matrix_width<3?m_pc_mass_matrix_width:2);
                for (int j(0); j<(dir+2*pc_mass_matrix_width); j++) {
                  for (int i(0); i<((1-dir)+2*pc_mass_matrix_width); i++) {
                    IntVect iL(ic);
                    iL[0] += (i-pc_mass_matrix_width);
                    iL[1] += (j-pc_mass_matrix_width);
                    int pL = (int) pmap_E_rhs[dit][dir](iL,0);
                    if (pL >= 0) {
                      auto val = sig_norm_nc * (dir == 0 ? m_sigma_zx_pc[dit](ic,idx)
                                                         : m_sigma_zy_pc[dit](ic,idx) );
                      insertOrAdd<int,Real>(this_row, pL, val);
                    }
                    idx++;
                  }
                }
              }
            }
#endif
          }


          CH_assert(this_row.size() <= m_pcmat_nnz);
          CH_assert(this_row.size() <= a_P.getNBands());

          std::vector<int> icols(0);
          std::vector<Real> vals(0);
          for (auto it = this_row.cbegin(); it != this_row.cend(); ++it) {
              icols.push_back( it->first );
              vals.push_back( it->second );
          }
          a_P.setRowValues( pc, this_row.size(), icols.data(), vals.data() );
        }
      }
    }
#endif

  }

  return;
}

#include "NamespaceFooter.H"

