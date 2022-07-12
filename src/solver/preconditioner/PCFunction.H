#ifndef _PCFUNCTION_H_
#define _PCFUNCTION_H_

#include "BandedMatrix.H"
#include "LinearFunction.H"

#include "NamespaceHeader.H"

template <class T, class Ops>
class PCFunction : public LinearFunction<T,Ops>
{
  public:
  
    virtual void update( const T& ) = 0;

    virtual void curTime(const double) = 0;
    virtual void curTimeStep(const double) = 0;

    virtual void setBaseSolution(const T&) = 0;
    virtual void setBaseRHS(const T&) = 0;

    virtual void setJFNKEps(Real) = 0;
    virtual void setIsLinear(bool) = 0;

    virtual const BandedMatrix& getBandedMatrixForm() const = 0;
    virtual BandedMatrix& getBandedMatrixForm() = 0;

    virtual bool isDefined() const = 0;

    virtual void setOutputIndent( const std::string& ) { }

  protected:

  private:
    
};

#include "NamespaceFooter.H"
#endif
