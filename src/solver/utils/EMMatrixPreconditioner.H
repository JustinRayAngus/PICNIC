#ifndef _EMMatrixPC_H_
#define _EMMatrixPC_H_

#include <string>

#include "Preconditioner.H"
#include "EMPCFunction.H"
#include "LinearFunctionSolverLibrary.H"

#define _EM_MATRIX_PC_ "matrixPC"

#include "NamespaceHeader.H"

template <class T, class Ops>
class EMMatrixPreconditioner : public Preconditioner<T,Ops>
{
  public:

    EMMatrixPreconditioner<T,Ops>()
    {
      m_is_defined = false;
      m_pc_mat = nullptr;
      m_pc_solver = nullptr;
      m_pc_solver_type = _LINSOLVER_GMRES_;
      m_op_indent = "";
      return;
    }

    virtual ~EMMatrixPreconditioner<T,Ops>()
    {
      delete m_pc_mat;
      delete m_pc_solver;
      return;
    }

    inline bool isDefined() const 
    { 
      return m_is_defined; 
    }

    virtual void define(  const T&, 
                          Ops* const,
                          const double );

    virtual void update( const T& a_X )
    {
      CH_TIME("EMMatrixPreconditioner::update()");
      if (!procID()) {
        // JRA: update called every time step with new mass matrix 
        // formalism ==> prints too much.
        //printf( "%s(%s) Updating preconditioner matrix.\n", 
        //        m_op_indent.c_str(), _EM_MATRIX_PC_ );
      }
      m_pc_mat->update(a_X);
      return;
    }

    virtual void apply(T& a_Y, const T& a_X)
    {
      CH_assert(isDefined());
      m_pc_solver->solve(a_Y, a_X);
      return;
    }

    inline void curTime(const double a_t) 
    { 
      m_pc_mat->curTime(a_t); 
    }

    inline void curTimeStep(const double a_dt) 
    { 
      m_pc_mat->curTimeStep(a_dt); 
    }

    virtual void setBaseSolution(const T& a_x) 
    { 
      m_pc_mat->setBaseSolution(a_x); 
    }

    virtual void setBaseRHS(const T& a_F) 
    { 
      m_pc_mat->setBaseRHS(a_F); 
    }

    inline void setJFNKEps(double a_eps) 
    {
      m_pc_mat->setJFNKEps(a_eps); 
    }

    inline void setIsLinear(bool a_b) 
    { 
      m_pc_mat->setIsLinear(a_b); 
    }

    virtual BandedMatrix& getBandedMatrix() 
    { 
      return(m_pc_mat->getBandedMatrixForm()); 
    }

    virtual void setOutputIndent( const std::string& a_indent ) 
    { 
      m_op_indent = a_indent;
      if (m_pc_mat) {
        m_pc_mat->setOutputIndent(a_indent);
      }
      if (m_pc_solver) {
        m_pc_solver->setOutputIndent(a_indent);
      }
    }

#ifdef with_petsc
    virtual void copyMatrixToPetsc( Mat a_P )
    {
      BandedMatrix& Pmat( m_pc_mat->getBandedMatrixForm() );
      Pmat.copyToPetscMat( a_P );
    }
#endif

  protected:
    
    bool m_is_defined;

    std::string m_pc_solver_type;
    std::string m_op_indent;

    Ops* m_ops;

    PCFunction<T,Ops>* m_pc_mat;
    LinearFunctionSolver<T,Ops>* m_pc_solver;

  private:

};

template <class T, class Ops>
void EMMatrixPreconditioner<T,Ops>::define( const T&      a_x,
                                            Ops* const    a_ops,
                                            const double  a_dtfac )
{
  CH_assert(!isDefined());

  m_ops = a_ops;

  m_pc_mat = new EMPCFunction<T,Ops>;
  m_pc_mat->define( a_x, m_ops, a_dtfac );

  ParmParse pp_pc(_EM_MATRIX_PC_);
  pp_pc.query("solver_type", m_pc_solver_type);

  if (m_pc_solver_type == _LINSOLVER_GMRES_ ) {
    m_pc_solver = new GMRESSolver<T,Ops>;
  } else if (m_pc_solver_type == _LINSOLVER_JACOBI_) {
    m_pc_solver = new PointJacobiSolver<T,Ops>;
  } else if (m_pc_solver_type == _LINSOLVER_GAUSSSEIDEL_) {
    m_pc_solver = new GaussSeidelSolver<T,Ops>;
  } else {
      MayDay::Error("Invalid choice of pc solver type");
  }

  m_pc_solver->define(m_pc_mat, a_x);
  m_is_defined = true;
  return;
}

#include "NamespaceFooter.H"
#endif

