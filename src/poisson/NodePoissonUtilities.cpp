#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "CONSTANTS.H"
#include "parstream.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "RealVect.H"
#include "NodePoissonUtilities.H"
#include "functionsF_F.H"
#include "PoissProbF_F.H"
#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "BoxIterator.H"
#include "BiCGStabSolver.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "NodeBCFunc.H"
#include "AMRNodeOp.H"


std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

void
locGetAcoef(LevelData<FArrayBox>& a_acoef,
            const DisjointBoxLayout& a_grids,
            const Real& a_dx)
{
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& pt = bit();
          RealVect loc(pt);
          loc += 0.5*RealVect::Unit;
          loc *= a_dx;
          Real val = 0;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              Real x = loc[idir];
              val += sin(x)*sin(x);
            }
        }
    }
}

/***************/
RealVect& getTrigRV()
{
  Real pi = 4.*atan(1.0);
  if (!GlobalBCRS::s_trigParsed)
    {
      GlobalBCRS::s_trigParsed = true;
      ParmParse pp;
      std::vector<Real> trigvec(SpaceDim);
      pp.getarr("trig",trigvec,0,SpaceDim);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          GlobalBCRS::s_trigvec[idir] = pi*trigvec[idir];
        }
    }
  return GlobalBCRS::s_trigvec;
}
void TrigValueNeum(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  RealVect gradPhi;
  FORT_GETGRADPHIPOINT(CHF_REALVECT(gradPhi),
                       CHF_CONST_REALVECT(trig),
                       CHF_CONST_REALVECT(xval));

  a_values[0] = gradPhi[*dir];
}

void ResistDiri(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  ParmParse pp;
  int whichMag;
  pp.get("which_mag", whichMag);
  for (int icomp = 0; icomp < 3; icomp++)
    {
      Real value;
      FORT_GETMAGPOINTRESIST(CHF_REAL(value),
                             CHF_CONST_REALVECT(trig),
                             CHF_CONST_REALVECT(xval),
                             CHF_INT(icomp),
                             CHF_INT(whichMag));
      a_values[icomp] = value;
    }
}

void ViscousDiri(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  ParmParse pp;
  int whichMag;
  pp.get("which_mag", whichMag);
  for (int icomp = 0; icomp < SpaceDim; icomp++)
    {
      Real value;
      FORT_GETMAGPOINTRESIST(CHF_REAL(value),
                             CHF_CONST_REALVECT(trig),
                             CHF_CONST_REALVECT(xval),
                             CHF_INT(icomp),
                             CHF_INT(whichMag));
      a_values[icomp] = value;
    }
}

void TrigValueDiri(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  Real value;
  FORT_GETPHIPOINT(CHF_REAL(value),
                   CHF_CONST_REALVECT(trig),
                   CHF_CONST_REALVECT(xval));
  a_values[0] = value;
}

void
nodeOutputData(const Vector<LevelData<NodeFArrayBox>* >& vectPhi,
               const Vector<DisjointBoxLayout>& vectGrids,
               const ProblemDomain&     domain,
               const Vector<int>& vectRatio,
               Real dxCoarsest,
               int numlevels,
               string filename,
               string varname)
{
  // for now, make a placeholders for these -- if more than
  // one comp, probably want to change vectNames accordingly
  Vector<string> vectNames(vectPhi[0]->nComp(), varname);
  Real dt = 1.0;
  Real time = 0.0;

#ifdef CH_USE_HDF5
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhi,
                        vectNames,
                        domain.domainBox(),
                        dxCoarsest,
                        dt, time,
                        vectRatio,
                        numlevels);
#endif

}

void
outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
           const Vector<DisjointBoxLayout>& vectGrids,
           const ProblemDomain&     domain,
           const Vector<int>& vectRatio,
           Real dxCoarsest,
           int numlevels,
           string filename,
           string varname)
{
#ifdef CH_USE_HDF5
  Vector<string> vectName(vectPhi[0]->nComp(), varname);
  Real time = 1;  //placeholder
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhi,
                        vectName,
                        domain.domainBox(),
                        dxCoarsest, time, time,
                        vectRatio,
                        numlevels);
#endif
}


void
outputVectorData(const Vector<LevelData<FArrayBox>* >& vectPhi,
                 const Vector<DisjointBoxLayout>& vectGrids,
                 const ProblemDomain&     domain,
                 const Vector<int>& vectRatio,
                 Real dxCoarsest,
                 int numlevels,
                 string filename,
                 Vector<string>& a_vectName)
{
#ifdef CH_USE_HDF5
  Real time = 1;  //placeholder
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhi,
                        a_vectName,
                        domain.domainBox(),
                        dxCoarsest, time, time,
                        vectRatio,
                        numlevels);
#endif
}

/*
  Set grid hierarchy from input file
*/
void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         PoissonParameters&         a_params)
{

  vectDomain.resize(a_params.numLevels);
  vectDx.resize(    a_params.numLevels);
  vectDx[0] = a_params.coarsestDx;
  for (int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
      vectDx[ilev] = vectDx[ilev-1]/a_params.refRatio[ilev-1];
    }


  vectDomain[0] = a_params.coarsestDomain;
  for (int ilev = 1;ilev < a_params.numLevels; ilev++)
    {
      vectDomain[ilev] = refine(vectDomain[ilev-1],a_params.refRatio[ilev-1]);
    }
}

int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             PoissonParameters&         a_params)
{
  Vector<ProblemDomain>     vectDomain;
  Vector<Real>              vectDx;
  getDomainsAndDxes(vectDomain, vectDx, a_params);

  int numlevels = a_params.numLevels;

  ParmParse pp("poisson");
  bool useEBGrids;
  // grid generation parameters

  vectGrids.resize(numlevels);
  bool readInGrids = false;
  pp.query("use_eb_grids", useEBGrids);
  if(pp.contains("read_in_grids")) pp.get("read_in_grids", readInGrids);

  if (readInGrids)
    {

      ProblemDomain levDomain = a_params.coarsestDomain;
      for (int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          Vector<Box>   boxes;
          char boxCountVar[100];
          int boxCount;
          sprintf(boxCountVar, "level_%d_box_count", ilev);
          pp.get(boxCountVar, boxCount);
          boxes.resize(boxCount);
          for (int ibox = 0; ibox < boxCount; ibox++)
            {
              char boxLoVar[100];
              char boxHiVar[100];
              sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
              sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
              Vector<int> boxLo, boxHi;
              pp.getarr(boxLoVar, boxLo, 0, SpaceDim);
              pp.getarr(boxHiVar, boxHi, 0, SpaceDim);
              IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));
              IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));
              boxes[ibox] = Box(ivLo, ivHi);
              if (!levDomain.contains(boxes[ibox]))
                {
                  MayDay::Error("box outside of domain");
                }
            }
          //check to see if level 0 domain is covered
          if (ilev == 0)
            {
              IntVectSet ivDom(levDomain.domainBox());
              for (int ibox = 0; ibox < boxes.size(); ibox++)
                {
                  ivDom -= boxes[ibox];
                }
              if (!ivDom.isEmpty())
                {
                  MayDay::Error("level 0 boxes must cover the domain");
                }
            }
          Vector<int>  proc(a_params.numLevels);
          LoadBalance(proc,boxes);
          vectGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
          levDomain.refine(a_params.refRatio[ilev]);
        }

    }
  else if (useEBGrids)
    {
      pout() << "all regular geometry" << endl;
      pout() << "ignoring grid parameters and making simple grids" << endl;

      Vector<int> proc(1, 0);
      Box coarBox = vectDomain[0].domainBox();
      Vector<Box> coarBoxes(1, coarBox);
      vectGrids[0] = DisjointBoxLayout(coarBoxes, proc, vectDomain[0]);

    }
  else
    {
      pout() << "using max_grid_size (" << a_params.maxGridSize << ") parameter to set boxes" << endl;
      
      Vector<Vector<Box> > oldBoxes(numlevels);
      Vector< Vector<int> > procAssign(numlevels);
      domainSplit(vectDomain[0], oldBoxes[0], a_params.maxGridSize, a_params.blockFactor);
      procAssign[0].resize(oldBoxes[0].size());
      LoadBalance(procAssign[0],oldBoxes[0]);

      vectGrids[0].define(oldBoxes[0],procAssign[0],vectDomain[0]);

    }
  return 0;
}




void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  ParmParse pp("poisson");
  Real bcVal;
  pp.get("bc_value",bcVal);
  a_values[0]=bcVal;
}


void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      if (!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp("poisson");
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          // don't do anything if periodic
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if (GlobalBCRS::s_bcLo[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const neum bcs lo for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig neum bcs lo for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueNeum,
                             i,
                             Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueDiri,
                             i,
                             Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 4)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "periodic bcs lo for direction " << i << endl;
                        }

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 5)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "reflective slip bcs lo for direction " << i << endl;
                        }
                      ReflectiveVectorBC(a_state,
                                         valid,
                                         a_dx,
                                         i,
                                         Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 6)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "no slip bcs lo for direction " << i << endl;
                        }
                      NoSlipVectorBC(a_state,
                                     valid,
                                     a_dx,
                                     i,
                                     Side::Lo, 1);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 7)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != 3)
                            {
                              MayDay::Error("bc function hardwired for ncomp == 3");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "Resistive diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ResistDiri,
                             i,
                             Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 8)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          if (a_state.nComp() != SpaceDim)
                            {
                              MayDay::Error("bc function hardwired for ncomp == SpaceDim");
                            }
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "Viscous diri bcs lo for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ViscousDiri,
                             i,
                             Side::Lo);

                    }
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                  if (GlobalBCRS::s_bcHi[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const neum bcs hi for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig neum bcs hi for direction " << i << endl;
                        }
                      NeumBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueNeum,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ParseValue,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 1)
                            {
                              MayDay::Error("using scalar bc function for vector");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             TrigValueDiri,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 4)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "periodic bcs hi for direction " << i << endl;
                        }
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 5)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "reflective slip bcs hi for direction " << i << endl;
                        }
                      ReflectiveVectorBC(a_state,
                                         valid,
                                         a_dx,
                                         i,
                                         Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 6)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() == 1)
                            {
                              MayDay::Error("using vector bc function for scalar");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "no slip bcs hi for direction " << i << endl;
                        }
                      NoSlipVectorBC(a_state,
                                     valid,
                                     a_dx,
                                     i,
                                     Side::Hi, 1);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 7)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != 3)
                            {
                              MayDay::Error("this bc hardwired for  ncomp=3");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "resisitive mhd diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ResistDiri,
                             i,
                             Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 8)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          if (a_state.nComp() != SpaceDim)
                            {
                              MayDay::Error("this bc hardwired for spacedim ncomp");
                            }
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "resisitive mhd diri bcs hi for direction " << i << endl;
                        }
                      DiriBC(a_state,
                             valid,
                             a_dx,
                             a_homogeneous,
                             ViscousDiri,
                             i,
                             Side::Hi);
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag hi");
                    }
                }
            } // end if is not periodic in ith direction
        }
    }
}

void NodeParseBC(NodeFArrayBox& a_state,
                 const Box& a_valid,
                 const ProblemDomain& a_domain,
                 Real a_dx,
                 bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      if (!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp("poisson");
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      Box valid = a_valid;
      for (int i=0; i<CH_SPACEDIM; ++i)
        {
          if (!a_domain.isPeriodic(i))
            {
              Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
              Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
              if (!a_domain.domainBox().contains(ghostBoxLo))
                {
                  if (GlobalBCRS::s_bcLo[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const neum bcs lo for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig neum bcs lo for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueNeum,
                                 i,
                                 Side::Lo);
                    }
                  else if (GlobalBCRS::s_bcLo[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "const diri bcs lo for direction "
                                 << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Lo);

                    }
                  else if (GlobalBCRS::s_bcLo[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatLo[i])
                        {
                          GlobalBCRS::s_printedThatLo[i] = true;
                          pout() << "trig diri bcs lo for direction "
                                 << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueDiri,
                                 i,
                                 Side::Lo);

                    }
                  else
                    {
                      MayDay::Error("bogus bc flag lo");
                    }
                }

              if (!a_domain.domainBox().contains(ghostBoxHi))
                {
                  if (GlobalBCRS::s_bcHi[i] == 1)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const neum bcs hi for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Hi);
                    }
                  if (GlobalBCRS::s_bcHi[i] == 2)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig neum bcs hi for direction "
                                 << i << endl;
                        }
                      NodeNeumBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueNeum,
                                 i,
                                 Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 0)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "const diri bcs hi for direction "
                                 << i << endl;
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 ParseValue,
                                 i,
                                 Side::Hi);
                    }
                  else if (GlobalBCRS::s_bcHi[i] == 3)
                    {
                      if (!GlobalBCRS::s_printedThatHi[i])
                        {
                          GlobalBCRS::s_printedThatHi[i] = true;
                          pout() << "trig diri bcs hi for direction "
                                 << i << endl; //
                        }
                      NodeDiriBC(a_state,
                                 valid,
                                 a_dx,
                                 a_homogeneous,
                                 TrigValueDiri,
                                 i,
                                 Side::Hi);
                    }
                  else
                    {
                      MayDay::Error("bogus bc flag hi");
                    }
                }
            } // end if periodic in ith direction
        }

    }

}


void
nodeDefineSolver( AMRMultiGrid<LevelData<NodeFArrayBox>>&  a_solver,
            const Vector<DisjointBoxLayout>&               a_grids,
                  LinearSolver<LevelData<NodeFArrayBox>>&  a_bottomSolver,
            const PoissonParameters&                       a_params )
{
  ParmParse pp2("poisson");
  AMRNodeOpFactory opFactory;
  opFactory.define(a_params.coarsestDomain,
                   a_grids,
                   a_params.refRatio,
                   a_params.coarsestDx,
                   &NodeParseBC, a_params.alpha, a_params.beta);

  ProblemDomain coarsestDomain(a_params.coarsestDomain);
  a_solver.define(coarsestDomain, opFactory,  &a_bottomSolver, a_params.numLevels);

  int numSmooth, numMG, maxIter;
  Real eps, hang;
  pp2.get("num_smooth", numSmooth);
  pp2.get("num_mg",     numMG);
  pp2.get("max_iterations", maxIter);
  pp2.get("tolerance", eps);
  pp2.get("hang",      hang);
  Real normThresh = 1.0e-30;
  
  // see section 4.4.4 in ChomboDesign doc
  a_solver.setSolverParameters(numSmooth, numSmooth, numSmooth,
                               numMG, maxIter, eps, hang, normThresh);
  a_solver.m_verbosity = 6;

}

/******/
int iscript(int icomp, int inorm, int ncomp)
{
  return icomp + inorm*ncomp;
}

/********/
void getCoarseLayoutsFromFine(Vector<DisjointBoxLayout>&       a_gridsCoar,
                              const Vector<DisjointBoxLayout>& a_gridsFine,
                              const PoissonParameters&         a_paramsCoar)
{
  int nlevels = a_paramsCoar.numLevels;
  a_gridsCoar.resize(nlevels);
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      CH_assert(a_gridsFine[ilev].coarsenable(2));
      coarsen(a_gridsCoar[ilev], a_gridsFine[ilev], 2);
    }

}
/********/
void PoissonParameters::coarsen(int a_factor)
{
  coarsestDx *= a_factor;
  coarsestDomain.coarsen(a_factor);
}
/********/
void PoissonParameters::refine(int a_factor)
{
  coarsestDx /= a_factor;
  coarsestDomain.refine(a_factor);
}
/********/
void getPoissonParameters( PoissonParameters&  a_params,
                     const DomainGrid&         a_mesh )
{
   ParmParse pp("poisson");
  
   a_params.refineThresh = 0.85; // only used with more than 1 level
   a_params.fillRatio = 0.5;     // only used with more then 1 level
   a_params.bufferSize = 1;      // not sure where this is used
   a_params.maxLevel = 0;        // only level 0 for no AMR
   a_params.numLevels = a_params.maxLevel + 1;
   
   pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
   pp.get("max_grid_size",a_params.maxGridSize); // box size for default dbl setup
   pp.get("block_factor",a_params.blockFactor);  // box size must be coarsenable by this factor

   CH_assert(a_params.maxGridSize >= a_params.blockFactor);

   const ProblemDomain& domain(a_mesh.getDomain());
   const RealVect& Xmin(a_mesh.getXmin());
   const RealVect& Xmax(a_mesh.getXmax());
   const IntVect nCells = domain.size();

   for (int idir = 0; idir < SpaceDim; idir++) {
      a_params.nCells[idir] = nCells[idir];
   }

   a_params.alpha = 0.0;
   a_params.beta = 1.0;
   a_params.verbosity = 3;

   IntVect lo = IntVect::Zero;
   IntVect hi = a_params.nCells;
   hi -= IntVect::Unit;

   Box crseDomBox(lo,hi);
   ProblemDomain crseDom(crseDomBox);
   for (int dir=0; dir<SpaceDim; dir++) {
      crseDom.setPeriodic(dir, domain.isPeriodic(dir));
   }
   a_params.coarsestDomain = crseDom;

   // set domain length
   std::vector<Real> dLArray(SpaceDim);
   for (int idir = 0; idir < SpaceDim; idir++) {
      a_params.domainLength[idir] = Xmax[idir] - Xmin[idir];
   }

   //derived stuff
   a_params.coarsestDx = a_params.domainLength[0]/a_params.nCells[0];
   a_params.probLo = RealVect::Zero;
   a_params.probHi = RealVect::Zero;
   a_params.probHi += a_params.domainLength;

}

