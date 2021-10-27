#ifndef _NONLINEAR_SOLVER_H_
#define _NONLINEAR_SOLVER_H_

#include "NonlinearFunction.H"

#include "NamespaceHeader.H"

template<class Vec, class Ops>
class NonlinearSolver
{
  public:

    NonlinearSolver<Vec,Ops>() { }
    virtual ~NonlinearSolver<Vec,Ops>() { }

    virtual void define(  const Vec&, 
                          Ops* const, 
                          NonlinearFunction<Vec,Ops>* const,
                          const double ) = 0;

    virtual void solve( Vec&,
                        const Vec&,
                        const double,
                        const double ) const = 0;

    virtual bool isDefined() const = 0;
    virtual void setParams(const double, const double, const int) = 0;

  protected:

  private:
};

#include "NamespaceFooter.H"
#endif
