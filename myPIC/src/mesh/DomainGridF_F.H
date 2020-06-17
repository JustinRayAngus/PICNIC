#ifndef _DOMAINGRIDF_F_H_
#define _DOMAINGRIDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGET_NC_MAPPED_COORDS 
#define GUARDGET_NC_MAPPED_COORDS 
// Prototype for Fortran procedure get_nc_mapped_coords ...
//
void FORTRAN_NAME( GET_NC_MAPPED_COORDS ,get_nc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) );

#define FORT_GET_NC_MAPPED_COORDS FORTRAN_NAME( inlineGET_NC_MAPPED_COORDS, inlineGET_NC_MAPPED_COORDS)
#define FORTNT_GET_NC_MAPPED_COORDS FORTRAN_NAME( GET_NC_MAPPED_COORDS, get_nc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_NC_MAPPED_COORDS, inlineGET_NC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_NC_MAPPED_COORDS");
 FORTRAN_NAME( GET_NC_MAPPED_COORDS ,get_nc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_NC_MAPPED_COORDS 

#ifndef GUARDGET_CC_MAPPED_COORDS 
#define GUARDGET_CC_MAPPED_COORDS 
// Prototype for Fortran procedure get_cc_mapped_coords ...
//
void FORTRAN_NAME( GET_CC_MAPPED_COORDS ,get_cc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) );

#define FORT_GET_CC_MAPPED_COORDS FORTRAN_NAME( inlineGET_CC_MAPPED_COORDS, inlineGET_CC_MAPPED_COORDS)
#define FORTNT_GET_CC_MAPPED_COORDS FORTRAN_NAME( GET_CC_MAPPED_COORDS, get_cc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_CC_MAPPED_COORDS, inlineGET_CC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_CC_MAPPED_COORDS");
 FORTRAN_NAME( GET_CC_MAPPED_COORDS ,get_cc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_CC_MAPPED_COORDS 

#ifndef GUARDGET_FC_MAPPED_COORDS 
#define GUARDGET_FC_MAPPED_COORDS 
// Prototype for Fortran procedure get_fc_mapped_coords ...
//
void FORTRAN_NAME( GET_FC_MAPPED_COORDS ,get_fc_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) );

#define FORT_GET_FC_MAPPED_COORDS FORTRAN_NAME( inlineGET_FC_MAPPED_COORDS, inlineGET_FC_MAPPED_COORDS)
#define FORTNT_GET_FC_MAPPED_COORDS FORTRAN_NAME( GET_FC_MAPPED_COORDS, get_fc_mapped_coords)

inline void FORTRAN_NAME(inlineGET_FC_MAPPED_COORDS, inlineGET_FC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_FC_MAPPED_COORDS");
 FORTRAN_NAME( GET_FC_MAPPED_COORDS ,get_fc_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_FC_MAPPED_COORDS 

#ifndef GUARDGET_EC_MAPPED_COORDS 
#define GUARDGET_EC_MAPPED_COORDS 
// Prototype for Fortran procedure get_ec_mapped_coords ...
//
void FORTRAN_NAME( GET_EC_MAPPED_COORDS ,get_ec_mapped_coords )(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) );

#define FORT_GET_EC_MAPPED_COORDS FORTRAN_NAME( inlineGET_EC_MAPPED_COORDS, inlineGET_EC_MAPPED_COORDS)
#define FORTNT_GET_EC_MAPPED_COORDS FORTRAN_NAME( GET_EC_MAPPED_COORDS, get_ec_mapped_coords)

inline void FORTRAN_NAME(inlineGET_EC_MAPPED_COORDS, inlineGET_EC_MAPPED_COORDS)(
      CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_FRA(xi) )
{
 CH_TIMELEAF("FORT_GET_EC_MAPPED_COORDS");
 FORTRAN_NAME( GET_EC_MAPPED_COORDS ,get_ec_mapped_coords )(
      CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_FRA(xi) );
}
#endif  // GUARDGET_EC_MAPPED_COORDS 

}

#endif
