#ifndef _Preconditioner_H_
#define _Preconditioner_H_

#include <string>
#ifdef with_petsc
#include <petscmat.h>
#endif

#include "NamespaceHeader.H"

template <class T, class Ops>
class Preconditioner
{
  public:

    Preconditioner() { }
    virtual ~Preconditioner() {;}

    virtual void define(const T&, 
                        Ops* const,
                        const double ) = 0;

    virtual void update( const T& ) = 0;

    virtual void apply(T&, const T&) = 0;
    virtual bool isDefined() const = 0;

    virtual void curTime(const double) = 0;
    virtual void curTimeStep(const double) = 0;

    /* additional stuff for matrix-free types */
    virtual void setBaseSolution(const T&) {}
    virtual void setBaseRHS(const T&) {}
    virtual void setJFNKEps(double a_e) {}
    virtual void setIsLinear(bool a_b) {}

    virtual void setOutputIndent( const std::string& ) { }

#ifdef with_petsc
    virtual void copyMatrixToPetsc( Mat ) = 0;
#endif

  protected:

  private:

};

typedef std::vector<int> DOFList;

#include "NamespaceFooter.H"
#endif

