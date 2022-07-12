/* GMRES Solver: Copied from Chombo's GMRES solver with some modifications
 * (the primary one being two template parameters instead of one)
 */

#ifndef _GMRESSOLVER_H_
#define _GMRESSOLVER_H_

#undef RECOMPUTE_FINAL_RESIDUAL
#undef RECOMPUTE_RESIDUAL
#undef CHECK_ORTHOGONALITY

#include <float.h>
#include "LinearFunctionSolver.H"
#include "ParmParse.H"
#include "parstream.H"

#define _LINSOLVER_GMRES_ "gmres"

#include "NamespaceHeader.H"

template <class Vec, class Ops>
class GMRESSolver : public LinearFunctionSolver<Vec,Ops>
{
  public:
  
    GMRESSolver<Vec,Ops>();
    ~GMRESSolver<Vec,Ops>();
  
    void clearData();
    void setConvergenceMetrics(double a_metric, double a_tolerance);

    void define(LinearFunction<Vec,Ops>* const, const Vec&);
    void solve(Vec& a_phi, const Vec& a_rhs, int a_its = -1);
  
    int restartLen() const { return m_restrtLen; }
    void setRestartLen(int mm);
    int getConvergedIterations() const { return m_its; }
    std::string getConvergedReason();
    int getExitStatus() const { return m_exitStatus; }

    inline void verbose(const bool a_v) 
    {
      m_verbose = a_v;
    }

    virtual void setOutputIndent( const std::string& a_indent ) 
    { 
      m_op_indent = a_indent;
      if (m_linear_function) {
        m_linear_function->setOutputIndent(m_op_indent);
      }
    }
    
    virtual int getConvergenceParams( int& a_last_iter, int& a_total_iter ) const
    { 
      a_last_iter = m_its;
      a_total_iter = m_count;
      return m_exitStatus;
    }

  private:

    void allocate();
    void CycleGMRES( Vec &xx, const Vec &bb,
                     int &reason, int &itcount, double &rnorm0,
                     const bool avoidnorms = false);
  
    void ResidualGMRES( Vec &a_vv, const Vec &a_xx,
                        const Vec &a_bb, Vec &a_temp );
  
    void BuildGMRESSoln( double nrs[], Vec &a_xx, const int it,
                         const Vec vv_0[] );
  
    void UpdateGMRESHessenberg( const int it, bool hapend, double &res );
    void TwoUnmodifiedGramSchmidtOrthogonalization( const int it );
    void ApplyAB( Vec &a_dest, const Vec &a_xx, Vec &a_temp ) const;
    void parseParameters(ParmParse&);
  
    bool m_verbose;

    double m_atol, m_rtol, m_small;
    int m_exitStatus, m_normType, m_maxits, m_its, m_restrtLen, m_count;
  
    LinearFunction<Vec,Ops>* m_linear_function;
    std::string m_op_indent;
  
    double *m_data;
    double *m_hes, *m_hh, *m_d, *m_ee, *m_dd, *m_bt;
    Vec    *m_work_arr;
};

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::allocate()
{
  int max_k         = m_restrtLen;
  int hh            = (max_k + 2) * (max_k + 1);
  int hes           = (max_k + 1) * (max_k + 1);
  int rs            = (max_k + 2);
  int cc            = (max_k + 1);
  int bt            = (max_k + 1);
  int size          = (hh + hes + rs + 2*cc + bt + 1);

  m_data = new double[size];
  m_hh = m_data;        // hh
  m_hes = m_data + hh;  // hes_
  m_d = m_hes + hes;    // rs_
  m_ee = m_d + rs;      // cc_
  m_dd = m_ee + cc;     // ss_
  m_bt = m_dd + cc;     // buffer for residual recalcs

  m_work_arr = 0;
}

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::clearData()
{
  delete [] m_data;
}

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::setRestartLen(int mm)
{
  CH_assert(mm>0);
  clearData();
  m_restrtLen = mm;
  allocate();
}

template <class Vec,class Ops>
GMRESSolver<Vec,Ops>::GMRESSolver()
  :m_verbose(true),
   m_atol(1.0E-12),
   m_rtol(1.0E-4),
   m_small(1.0E-30),
   m_exitStatus(-1),
   m_normType(2),
   m_maxits(1000),
   m_restrtLen(30),
   m_count(0),
   m_linear_function(nullptr),
   m_op_indent("")
{
  allocate();
}

template <class Vec,class Ops>
GMRESSolver<Vec,Ops>::~GMRESSolver()
{
  clearData();
}

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::define(  LinearFunction<Vec,Ops>* const a_func,
                                    const Vec& )
{
  ParmParse ppGmres ( "gmres" );
  parseParameters( ppGmres  );

  m_linear_function = a_func;
  return;
}

template <class Vec, class Ops>
void GMRESSolver<Vec,Ops>::parseParameters( ParmParse& a_pp)
{
  a_pp.query("verbose", m_verbose); /* verbosity */
  a_pp.query("atol",    m_atol);    /* absolute tolerance */
  a_pp.query("rtol",    m_rtol);    /* relative tolerance */
  a_pp.query("maxits",  m_maxits);  /* maximum iterations */
}

///
/**         m_exitStatus
   set = -1 if solver exited for an unknown reason
   set =  0 if solver converged to tolerance
   set =  1 if max number of iterations exceeded
   set =  2 if rho = 0
 */
template <class Vec, class Ops>
std::string GMRESSolver<Vec,Ops>::getConvergedReason()
{
  std::string convergedReason;
  if      (m_exitStatus == -1) convergedReason = "unknown";
  else if (m_exitStatus ==  0) convergedReason = "norm less than absolute or relative tolerance";
  else if (m_exitStatus ==  1) convergedReason = "number of iterations exceeded maximum";
  else if (m_exitStatus ==  2) convergedReason = "rho = 0";
  return convergedReason;
}

#define HH(a,b)  (m_hh        + (b)*(m_restrtLen+2) + (a))
#define HES(a,b) (m_hes       + (b)*(m_restrtLen+1) + (a))
#define CC(a)    (m_ee        + (a))
#define SS(a)    (m_dd        + (a))
#define GRS(a)   (m_d         + (a))

/* vector names */
#define VEC_OFFSET 3
#define VEC_TEMP_RHS       m_work_arr[0]
#define VEC_TEMP_LHS       m_work_arr[1]
#define VEC_TEMP           m_work_arr[2]
#define VEC_VV(i)          m_work_arr[VEC_OFFSET + i]

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::solve( Vec& a_xx, const Vec& a_bb, int a_its )
{
  CH_TIME("GMRESSolver::solve()");

  if (a_its < 0) a_its = m_maxits;

  const int nwork = VEC_OFFSET + m_restrtLen + 1; // flex = VEC_OFFSET + 2*(m_restrtLen + 1);
  m_work_arr = new Vec[nwork];
  m_linear_function->create(VEC_TEMP_RHS, a_bb);
  m_linear_function->create(VEC_TEMP_LHS, a_xx);
  m_linear_function->create(VEC_TEMP,     a_xx);
  for (int i=VEC_OFFSET;i<nwork;i++)
    {
      m_linear_function->create(m_work_arr[i], a_bb);
    }

  double rnorm0 = 0.0;
  Vec &vv_0 = VEC_VV(0);

  /* Compute the initial (preconditioned) residual (into 'vv_0')*/
  m_linear_function->assign( vv_0, a_bb );
  m_linear_function->setToZero( a_xx );

  // doit
  m_its = 0; m_exitStatus = -1;
  CycleGMRES( a_xx, a_bb, m_exitStatus, m_its, rnorm0 );
  // loop for restarts
  while ( m_exitStatus==-1 && m_its < a_its )
  {
    ResidualGMRES( vv_0, a_xx, a_bb, VEC_TEMP_RHS );
    CycleGMRES( a_xx, a_bb, m_exitStatus, m_its, rnorm0 );
  }
  if (m_exitStatus==-1 && m_its >= a_its)  m_exitStatus = 1;

#if defined(RECOMPUTE_FINAL_RESIDUAL) && !defined(RECOMPUTE_RESIDUAL)
  {
     Vec &Mb = VEC_TEMP_LHS;
     ResidualGMRES(Mb, a_xx, a_bb, VEC_TEMP_RHS);

     double res_norm = m_linear_function->norm( Mb, m_normType );
     if (procID()==0) printf("            Recalculated final absolute residual = %1.6e\n", res_norm);
  }
#endif
  
  // clean up
  for (int i=0;i<nwork;i++)
    {
      m_linear_function->clear(m_work_arr[i]);
    }
  delete [] m_work_arr; m_work_arr = 0;

  m_count += m_its;
}

#define CONVERGED(r0,r) (r<r0*m_rtol || r<m_atol)

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::CycleGMRES( Vec &a_xx, const Vec &a_bb,
                 int &a_reason, int &a_itcount, double &a_rnorm0,
                                 const bool a_avoidnorms/*=false*/ )
{
  CH_TIME("GMRESSolver::CycleGMRES()");

  CH_assert(m_linear_function != 0);
  double        res,hapbnd,tt;
  int         it;
  bool        hapend = false;
  Vec           &vv_0 = VEC_VV(0);

  /* scale VEC_VV (the initial residual) */
  res = m_linear_function->norm( vv_0, m_normType );
  *GRS(0) = res;

  /* check for the convergence */
  if ( res == 0. )
  {
    a_reason = 0; // 0 == converged,
    return;
  }
  // normilize
  tt = 1./res;
  m_linear_function->scale( vv_0, tt);

  if ( a_itcount == 0 ) a_rnorm0 = res;
  bool conv = CONVERGED( a_rnorm0, res );
  a_reason = conv ? 0 : -1;
  it = 0;
  while ( it < m_restrtLen && a_itcount < m_maxits )
  {
#ifndef RECOMPUTE_RESIDUAL
    if ( (it!=0 || a_itcount==0) )
    {
      if (m_verbose && !procID()) {
        printf("%s(%s) iter  = %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
                m_op_indent.c_str(), _LINSOLVER_GMRES_, a_itcount, res, res / a_rnorm0);
      }
    }
#endif

    if (a_reason == 0) break;

    const Vec &vv_it = VEC_VV(it);
    Vec &vv_it1 = VEC_VV(it+1);

    // apply AB
    {
      Vec &Mb = VEC_TEMP_LHS;
      ApplyAB( vv_it1, vv_it, Mb );
    }
    /* update hessenberg matrix and do Gram-Schmidt */
    TwoUnmodifiedGramSchmidtOrthogonalization(it);

    /* vv(i+1) . vv(i+1) */
    tt = m_linear_function->norm( vv_it1, m_normType );
    /* check for the happy breakdown */
    hapbnd = 1.e-99; // hard wired happy tol!!!
    if (tt < hapbnd)
    {
      hapend = true;
    }
    else
    {
      m_linear_function->scale( vv_it1, 1./tt);
    }

#ifdef CHECK_ORTHOGONALITY
    for (int j=0; j<it+1; ++j) {
       const Vec &vv_j = VEC_VV(j);
       double prod;
       m_linear_function->mDotProduct(vv_it1, 1, &vv_j, &prod);
       if (procID()==0) cout << string(18,' ') << it+1 << " " << j << " dot product = " << prod << endl;
    }
#endif

    /* save the magnitude */
    *HH(it+1,it) = tt; *HES(it+1,it) = tt;

    UpdateGMRESHessenberg( it, hapend, res );

#ifdef RECOMPUTE_RESIDUAL
    {
       m_linear_function->setToZero(VEC_TEMP);

       BuildGMRESSoln( m_bt, VEC_TEMP, it, &VEC_VV(0) );

       Vec &Mb = VEC_TEMP_LHS;
       ResidualGMRES(Mb, VEC_TEMP, a_bb, VEC_TEMP_RHS);

       double recalc_norm = m_linear_function->norm( Mb, m_normType );
       if (procID()==0) printf("                  Absolute residual = %1.6e (alg), %1.6e (recalc)\n", res, recalc_norm);

       res = recalc_norm;
    }
#endif
    
    // increment
    it++; (a_itcount)++;
    conv = CONVERGED( a_rnorm0, res );
    a_reason = conv ? 0 : -1;
    /* Catch error in happy breakdown and signal convergence and break from loop */
    if ( hapend )
    {
      break;
    }
  }
  
  /* Monitor if we know that we will not return for a restart */
  if ( (a_reason!=0 || a_itcount>=m_maxits) )
  {
    if (m_verbose && !procID()) {
      printf("%s(%s) iter  = %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
              m_op_indent.c_str(), _LINSOLVER_GMRES_, a_itcount, res, res / a_rnorm0);
    }
  }

  /*
    Down here we have to solve for the "best" coefficients of the Krylov
    columns, add the solution values together, and possibly unwind the
    preconditioning from the solution
  */
  /* Form the solution (or the solution so far) */
  BuildGMRESSoln( GRS(0), a_xx, it-1, &VEC_VV(0) );
}

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::ResidualGMRES( Vec &a_vv, const Vec &a_xx,
                                    const Vec &a_bb, Vec &a_temp_rhs )
{
  CH_TIME("GMRESSolver::ResidualGMRES()");

  CH_assert(m_linear_function != 0);
  m_linear_function->applyOp( a_temp_rhs, a_xx, m_its); /* Vec <- Ax */
  m_linear_function->assign( a_vv, a_bb );
  m_linear_function->incr( a_vv, a_temp_rhs, -1.0); /* b - A(B)x */
}

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::BuildGMRESSoln( double nrs[], Vec &a_xx, const int it,
                                     const Vec vv_0[] )
{
  double tt;
  int ii,k,j;

  /* Solve for solution vector that minimizes the residual */

  /* If it is < 0, no gmres steps have been performed */
  if (it < 0)
  {
    return;
  }
  if (*HH(it,it) != 0.0)
  {
    nrs[it] = *GRS(it) / *HH(it,it);
  }
  else
  {
    nrs[it] = 0.0;
  }

  for (ii=1; ii<=it; ii++)
  {
    k   = it - ii;
    tt  = *GRS(k);
    for (j=k+1; j<=it; j++) tt  = tt - *HH(k,j) * nrs[j];
    nrs[k]   = tt / *HH(k,k);
  }

  /* Accumulate the correction to the solution of the preconditioned problem in TEMP */
  Vec &temp = VEC_TEMP_RHS;
  m_linear_function->setToZero(temp);
  //temp->MAXPY( it+1, nrs, vv_0 );
  for (ii=0; ii<it+1; ii++)
    {
      m_linear_function->incr(temp, vv_0[ii], nrs[ii]);
    }
  /* unwind pc */
  /*If we preconditioned on the right, we need to solve for the correction to
    the unpreconditioned problem */
  Vec &temp_matop = VEC_TEMP_LHS;
  //ierr = pc->Apply( temp, temp_matop );CHKERRQ(ierr);
  m_linear_function->preCond( temp_matop, temp );
  m_linear_function->incr( a_xx, temp_matop, 1.0 );
}

/* GMRESSolver::UpdateGMRESHessenberg *****************************************
 *
 *   INPUT:
 *     - it:
 *     - hapend: good breakdown?
 *     - res: residual (out)
 *
 *   SIDE EFFECTS:
 *     - sets 'nwork_' and allocs 'work_'
 *
 *   RETURN:
 *     - PETSc error code
 */
template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::UpdateGMRESHessenberg( const int it, bool hapend, double &res )
{
  CH_TIME("GMRESSolver::UpdateGMRESHessenberg()");
  double *hh,*cc,*ss,tt;
  int   j;

  hh  = HH(0,it);
  cc  = CC(0);
  ss  = SS(0);

  /* Apply all the previously computed plane rotations to the new column
     of the Hessenberg matrix */
  for (j=1; j<=it; j++)
  {
    tt  = *hh;
    *hh = *cc * tt + *ss * *(hh+1);
    hh++;
    *hh = *cc++ * *hh - (*ss++ * tt);
  }

  /*
    compute the new plane rotation, and apply it to:
     1) the right-hand-side of the Hessenberg system
     2) the new column of the Hessenberg matrix
     thus obtaining the updated value of the residual
  */
  if ( !hapend )
  {
    tt = sqrt( *hh * *hh + *(hh+1) * *(hh+1) );
    *cc       = *hh / tt;
    *ss       = *(hh+1) / tt;
    *GRS(it+1) = - (*ss * *GRS(it));
    *GRS(it)   = *cc * *GRS(it);
    *hh       = *cc * *hh + *ss * *(hh+1);
    res      = Abs( *GRS(it+1) );
  }
  else
  {
    /* happy breakdown: HH(it+1, it) = 0, therfore we don'Vec need to apply
       another rotation matrix (so RH doesn'Vec change).  The new residual is
       always the new sine term times the residual from last time (GRS(it)),
       but now the new sine rotation would be zero...so the residual should
       be zero...so we will multiply "zero" by the last residual.  This might
       not be exactly what we want to do here -could just return "zero". */

    res = 0.0;
  }
}

/*
  uses 1 iteration of iterative refinement of UNMODIFIED Gram-Schmidt.
  It can give better performance when running in a parallel
  environment and in some cases even in a sequential environment (because
  MAXPY has more data reuse).

  Care is taken to accumulate the updated HH/HES values.
 */
template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::TwoUnmodifiedGramSchmidtOrthogonalization( const int it )
{
  CH_TIME("GMRESSolver::TwoUnmodifedGramSchmidtOrthogonalization()");
  double     *hh,*hes,*lhh = 0;
  Vec        &vv_1 = VEC_VV(it+1);
  const Vec  *vv_0 = &(VEC_VV(0));

  /* Don'Vec allocate small arrays */
  lhh = new double[it+1];

  /* update Hessenberg matrix and do unmodified Gram-Schmidt */
  hh  = HH(0,it);
  hes = HES(0,it);

  /* Clear hh and hes since we will accumulate values into them */
  for (int j=0; j<=it; j++)
  {
    hh[j]  = 0.0;
    hes[j] = 0.0;
  }

  for ( int ncnt = 0 ; ncnt < 2 ; ncnt++ )
  {
    /*
       This is really a matrix-vector product, with the matrix stored
       as pointer to rows
    */
    m_linear_function->mDotProduct(vv_1, it+1, vv_0, lhh);

    /*
      This is really a matrix vector product:
      [h[0],h[1],...]*[ v[0]; v[1]; ...] subtracted from v[it+1].
    */
    for (int j=0; j<=it; j++) lhh[j] = - lhh[j];
    for (int j=0; j<=it; j++)
      {
        m_linear_function->incr(vv_1, vv_0[j], lhh[j]);
      }
    for (int j=0; j<=it; j++)
    {
      hh[j]  -= lhh[j];     /* hh += <v,vnew> */
      hes[j] += lhh[j];     /* hes += - <v,vnew> */
    }
  }

  delete [] lhh;
}

/* PromSolver::PromPCApplyBAorAB ******************************************
 *
 */
template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::ApplyAB( Vec &a_dest, const Vec &a_xx, Vec &a_tmp_lhs ) const
{
  CH_TIME("GMRESSolver::ApplyAB()");
  m_linear_function->preCond( a_tmp_lhs, a_xx );
  m_linear_function->applyOp( a_dest, a_tmp_lhs, m_its);
}

template <class Vec,class Ops>
void GMRESSolver<Vec,Ops>::setConvergenceMetrics(double a_metric,
                                              double a_tolerance)
{
  m_atol = a_tolerance;
}

#include "NamespaceFooter.H"
#endif /*_GMRESSOLVER_H_*/
