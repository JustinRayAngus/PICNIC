#ifndef _LINEARFUNCTIONSOLVER_H_
#define _LINEARFUNCTIONSOLVER_H_

#include <string>
#include "LinearFunction.H"

#include "NamespaceHeader.H"

template <class Vec, class Ops>
class LinearFunctionSolver
{
  public:

    virtual ~LinearFunctionSolver() {}

    virtual void define(  LinearFunction<Vec,Ops>* const, 
                          const Vec&) = 0;

    virtual void solve( Vec& a_x, 
                        const Vec& a_r, 
                        int a_its = -1) = 0;

    virtual void verbose( const bool ) = 0;

    virtual void setOutputIndent( const std::string& ) { }
    
    virtual int getConvergenceParams( int& a_last_iter, int& a_total_iter ) const
    { 
      a_last_iter = 0;
      a_total_iter = 0;
      return 0;
    }

};

#include "NamespaceFooter.H"
#endif
