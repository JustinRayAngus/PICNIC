#ifndef _NonlinearFunction_H_
#define _NonlinearFunction_H_

#include "LinearFunction.H"

#include "NamespaceHeader.H"

template <class Vec, class Ops>
class NonlinearFunction
{
  public:
    virtual ~NonlinearFunction() {}
    virtual void evalFunction(Vec&, const Vec&, const int) = 0;
    virtual void define(const Vec&, Ops* const, const double, const bool) = 0;
    virtual LinearFunction<Vec,Ops>& getJacobian() = 0;
    virtual void setRHS( const Vec& ) = 0;

    virtual void curTime( const double ) { }
    virtual void curTimeStep( const double ) { }

    virtual bool isDefined() const = 0;

  protected:

  private:

};

#include "NamespaceFooter.H"
#endif
