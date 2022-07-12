#ifndef _BandedMatrixGaussSeidelSolver_H_
#define _BandedMatrixGaussSeidelSolver_H_

#include "BandedMatrixSolver.H"

#define _BANDED_MATRIX_SOLVER_GAUSS_SEIDEL_ "gauss_seidel"

#include "NamespaceHeader.H"

template <class Vec>
class BandedMatrixGaussSeidelSolver : public BandedMatrixSolver<Vec>
{
  public:

    BandedMatrixGaussSeidelSolver<Vec>()
    { 
      m_is_defined = false; 
    }

    ~BandedMatrixGaussSeidelSolver<Vec>() {}

    void define( const Vec& );
    void solve( Vec&, const BandedMatrix&, const Vec&, int a_its = -1);

    bool isDefined() const { return m_is_defined; }

    inline void setOutputIndent( const std::string& a_indent )
    {
      m_op_indent = a_indent;
    }

  private:
    
    bool m_is_defined;
    bool m_verbose;
    int  m_its; /*!< Number of GaussSeidel iterations */
    std::unordered_map<int,double> m_nonlocal_data;
    std::string m_op_indent;

    Vec m_Z;
    void parseParameters(ParmParse&);
};

template <class Vec>
void BandedMatrixGaussSeidelSolver<Vec>::define( const Vec& a_U )
{
  /* default values */
  m_its = 100;
  m_verbose = false;

  ParmParse pp_GS( _BANDED_MATRIX_SOLVER_GAUSS_SEIDEL_ );
  parseParameters( pp_GS );

  if (m_verbose) m_Z.define( a_U );
  m_is_defined = true;
}

template <class Vec>
void BandedMatrixGaussSeidelSolver<Vec>::parseParameters( ParmParse& a_pp )
{
  a_pp.query("its", m_its);
  a_pp.query("verbose", m_verbose);
  return;
}

template <class Vec>
void BandedMatrixGaussSeidelSolver<Vec>::solve( Vec&                a_x, 
                                                const BandedMatrix& a_A, 
                                                const Vec&          a_b, 
                                                int                 a_its)
{
  CH_assert(m_is_defined);
  CH_assert(a_A.isAssembled());

  if (a_its < 0) a_its = m_its;

  for (int n(0); n < a_its; n++) {

    if (m_verbose) {
      BandedMatrixSolver<Vec>::applyOp(m_Z, a_x, a_A);
      double err_norm = (a_b-m_Z).norm();
      if (!procID()) {
        printf("%s(%s) iter = %3d, residual = %1.6e\n",
                m_op_indent.c_str(),
                _BANDED_MATRIX_SOLVER_GAUSS_SEIDEL_,
                n, err_norm );
      }
    }

    BandedMatrixSolver<Vec>::exchangeData(a_A, a_x, m_nonlocal_data);

    for (int i(0); i < a_A.nRows(); i++) {

      double *data_ptr      = a_A.data()  + i*a_A.getNBands();
      int  *col_ptr       = a_A.iCols() + i*a_A.getNBands();

      double diag;
      double rhs = a_b[i];

      for (int k(0); k < a_A.nCols(i); k++) {

        int icol = col_ptr[k];

        if (a_A.isLocal(icol)) {
          
          int j = a_A.globalToLocal(icol);
          if (j == i) {
            diag = data_ptr[k];
          } else {
            rhs -= ( data_ptr[k] * a_x[j] );
          }

        } else {

          rhs -= (data_ptr[k] * m_nonlocal_data.at(icol) );

        }
      }

      a_x[i] = rhs / diag;
    }
  }

  return;
}

#include "NamespaceFooter.H"
#endif
