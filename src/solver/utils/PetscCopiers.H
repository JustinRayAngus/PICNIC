#ifndef _PETSCCOPIERS_H_
#define _PETSCCOPIERS_H_

#ifdef with_petsc

#include "NamespaceHeader.H"

namespace PetscCopiers{

  template<class T>
  void cpAppVecToPetsc( Vec&     a_Y, 
                        const T& a_U )
  {
    PetscScalar *Yarr;
    VecGetArray(a_Y,&Yarr);
    a_U.copyTo( (Real* const) Yarr );
    VecRestoreArray(a_Y,&Yarr);
  }
  
  template<class T>
  void cpAppVecFromPetsc( const Vec&  a_Y, 
                          T&          a_U )
  {
    const PetscScalar *Yarr;
    VecGetArrayRead(a_Y,&Yarr);
    a_U.copyFrom( (const Real* const) Yarr );
    VecRestoreArrayRead(a_Y,&Yarr);
  }

}

#include "NamespaceFooter.H"
#endif
#endif
