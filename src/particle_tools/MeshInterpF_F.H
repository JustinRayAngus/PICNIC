#ifndef _MESHINTERPF_F_H_
#define _MESHINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMOMENT_DEPOSIT 
#define GUARDMOMENT_DEPOSIT 
// Prototype for Fortran procedure moment_deposit ...
//
void FORTRAN_NAME( MOMENT_DEPOSIT ,moment_deposit )(
      CHFp_FRA1(moment)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(kernal)
      ,CHFp_CONST_REAL(q) );

#define FORT_MOMENT_DEPOSIT FORTRAN_NAME( inlineMOMENT_DEPOSIT, inlineMOMENT_DEPOSIT)
#define FORTNT_MOMENT_DEPOSIT FORTRAN_NAME( MOMENT_DEPOSIT, moment_deposit)

inline void FORTRAN_NAME(inlineMOMENT_DEPOSIT, inlineMOMENT_DEPOSIT)(
      CHFp_FRA1(moment)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(kernal)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_MOMENT_DEPOSIT");
 FORTRAN_NAME( MOMENT_DEPOSIT ,moment_deposit )(
      CHFt_FRA1(moment)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(kernal)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDMOMENT_DEPOSIT 

#ifndef GUARDMOMENT3V_DEPOSIT 
#define GUARDMOMENT3V_DEPOSIT 
// Prototype for Fortran procedure moment3v_deposit ...
//
void FORTRAN_NAME( MOMENT3V_DEPOSIT ,moment3v_deposit )(
      CHFp_FRA(energy)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(kernal0)
      ,CHFp_CONST_REAL(kernal1)
      ,CHFp_CONST_REAL(kernal2)
      ,CHFp_CONST_REAL(w) );

#define FORT_MOMENT3V_DEPOSIT FORTRAN_NAME( inlineMOMENT3V_DEPOSIT, inlineMOMENT3V_DEPOSIT)
#define FORTNT_MOMENT3V_DEPOSIT FORTRAN_NAME( MOMENT3V_DEPOSIT, moment3v_deposit)

inline void FORTRAN_NAME(inlineMOMENT3V_DEPOSIT, inlineMOMENT3V_DEPOSIT)(
      CHFp_FRA(energy)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(kernal0)
      ,CHFp_CONST_REAL(kernal1)
      ,CHFp_CONST_REAL(kernal2)
      ,CHFp_CONST_REAL(w) )
{
 CH_TIMELEAF("FORT_MOMENT3V_DEPOSIT");
 FORTRAN_NAME( MOMENT3V_DEPOSIT ,moment3v_deposit )(
      CHFt_FRA(energy)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(kernal0)
      ,CHFt_CONST_REAL(kernal1)
      ,CHFt_CONST_REAL(kernal2)
      ,CHFt_CONST_REAL(w) );
}
#endif  // GUARDMOMENT3V_DEPOSIT 

#ifndef GUARDSET_PARTICLE_WEIGHT 
#define GUARDSET_PARTICLE_WEIGHT 
// Prototype for Fortran procedure set_particle_weight ...
//
void FORTRAN_NAME( SET_PARTICLE_WEIGHT ,set_particle_weight )(
      CHFp_REAL(weight)
      ,CHFp_CONST_REAL(partsPerCell)
      ,CHFp_CONST_FRA1(density)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_SET_PARTICLE_WEIGHT FORTRAN_NAME( inlineSET_PARTICLE_WEIGHT, inlineSET_PARTICLE_WEIGHT)
#define FORTNT_SET_PARTICLE_WEIGHT FORTRAN_NAME( SET_PARTICLE_WEIGHT, set_particle_weight)

inline void FORTRAN_NAME(inlineSET_PARTICLE_WEIGHT, inlineSET_PARTICLE_WEIGHT)(
      CHFp_REAL(weight)
      ,CHFp_CONST_REAL(partsPerCell)
      ,CHFp_CONST_FRA1(density)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_SET_PARTICLE_WEIGHT");
 FORTRAN_NAME( SET_PARTICLE_WEIGHT ,set_particle_weight )(
      CHFt_REAL(weight)
      ,CHFt_CONST_REAL(partsPerCell)
      ,CHFt_CONST_FRA1(density)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDSET_PARTICLE_WEIGHT 

#ifndef GUARDNGP_DEPOSIT 
#define GUARDNGP_DEPOSIT 
// Prototype for Fortran procedure ngp_deposit ...
//
void FORTRAN_NAME( NGP_DEPOSIT ,ngp_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_NGP_DEPOSIT FORTRAN_NAME( inlineNGP_DEPOSIT, inlineNGP_DEPOSIT)
#define FORTNT_NGP_DEPOSIT FORTRAN_NAME( NGP_DEPOSIT, ngp_deposit)

inline void FORTRAN_NAME(inlineNGP_DEPOSIT, inlineNGP_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_NGP_DEPOSIT");
 FORTRAN_NAME( NGP_DEPOSIT ,ngp_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDNGP_DEPOSIT 

#ifndef GUARDNGP_INTERPOLATE 
#define GUARDNGP_INTERPOLATE 
// Prototype for Fortran procedure ngp_interpolate ...
//
void FORTRAN_NAME( NGP_INTERPOLATE ,ngp_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_NGP_INTERPOLATE FORTRAN_NAME( inlineNGP_INTERPOLATE, inlineNGP_INTERPOLATE)
#define FORTNT_NGP_INTERPOLATE FORTRAN_NAME( NGP_INTERPOLATE, ngp_interpolate)

inline void FORTRAN_NAME(inlineNGP_INTERPOLATE, inlineNGP_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_NGP_INTERPOLATE");
 FORTRAN_NAME( NGP_INTERPOLATE ,ngp_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDNGP_INTERPOLATE 

#ifndef GUARDCIC_DEPOSIT 
#define GUARDCIC_DEPOSIT 
// Prototype for Fortran procedure cic_deposit ...
//
void FORTRAN_NAME( CIC_DEPOSIT ,cic_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_INTVECT(stag)
      ,CHFp_CONST_REAL(kernal)
      ,CHFp_CONST_REAL(q) );

#define FORT_CIC_DEPOSIT FORTRAN_NAME( inlineCIC_DEPOSIT, inlineCIC_DEPOSIT)
#define FORTNT_CIC_DEPOSIT FORTRAN_NAME( CIC_DEPOSIT, cic_deposit)

inline void FORTRAN_NAME(inlineCIC_DEPOSIT, inlineCIC_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_INTVECT(stag)
      ,CHFp_CONST_REAL(kernal)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_CIC_DEPOSIT");
 FORTRAN_NAME( CIC_DEPOSIT ,cic_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp)
      ,CHFt_CONST_INTVECT(stag)
      ,CHFt_CONST_REAL(kernal)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDCIC_DEPOSIT 

#ifndef GUARDCIC_INTERPOLATE 
#define GUARDCIC_INTERPOLATE 
// Prototype for Fortran procedure cic_interpolate ...
//
void FORTRAN_NAME( CIC_INTERPOLATE ,cic_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_CIC_INTERPOLATE FORTRAN_NAME( inlineCIC_INTERPOLATE, inlineCIC_INTERPOLATE)
#define FORTNT_CIC_INTERPOLATE FORTRAN_NAME( CIC_INTERPOLATE, cic_interpolate)

inline void FORTRAN_NAME(inlineCIC_INTERPOLATE, inlineCIC_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_CIC_INTERPOLATE");
 FORTRAN_NAME( CIC_INTERPOLATE ,cic_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDCIC_INTERPOLATE 

#ifndef GUARDCIC_DEPOSIT_CURRENT 
#define GUARDCIC_DEPOSIT_CURRENT 
// Prototype for Fortran procedure cic_deposit_current ...
//
void FORTRAN_NAME( CIC_DEPOSIT_CURRENT ,cic_deposit_current )(
      CHFp_FRA1(rhoVx)
      ,CHFp_FRA1(rhoVy)
      ,CHFp_FRA1(rhoVz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_REAL(vpx)
      ,CHFp_CONST_REAL(vpy)
      ,CHFp_CONST_REAL(vpz)
      ,CHFp_CONST_REAL(qp) );

#define FORT_CIC_DEPOSIT_CURRENT FORTRAN_NAME( inlineCIC_DEPOSIT_CURRENT, inlineCIC_DEPOSIT_CURRENT)
#define FORTNT_CIC_DEPOSIT_CURRENT FORTRAN_NAME( CIC_DEPOSIT_CURRENT, cic_deposit_current)

inline void FORTRAN_NAME(inlineCIC_DEPOSIT_CURRENT, inlineCIC_DEPOSIT_CURRENT)(
      CHFp_FRA1(rhoVx)
      ,CHFp_FRA1(rhoVy)
      ,CHFp_FRA1(rhoVz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_REAL(vpx)
      ,CHFp_CONST_REAL(vpy)
      ,CHFp_CONST_REAL(vpz)
      ,CHFp_CONST_REAL(qp) )
{
 CH_TIMELEAF("FORT_CIC_DEPOSIT_CURRENT");
 FORTRAN_NAME( CIC_DEPOSIT_CURRENT ,cic_deposit_current )(
      CHFt_FRA1(rhoVx)
      ,CHFt_FRA1(rhoVy)
      ,CHFt_FRA1(rhoVz)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp)
      ,CHFt_CONST_REAL(vpx)
      ,CHFt_CONST_REAL(vpy)
      ,CHFt_CONST_REAL(vpz)
      ,CHFt_CONST_REAL(qp) );
}
#endif  // GUARDCIC_DEPOSIT_CURRENT 

#ifndef GUARDTSC_DEPOSIT_CURRENT 
#define GUARDTSC_DEPOSIT_CURRENT 
// Prototype for Fortran procedure tsc_deposit_current ...
//
void FORTRAN_NAME( TSC_DEPOSIT_CURRENT ,tsc_deposit_current )(
      CHFp_FRA1(rhoVx)
      ,CHFp_FRA1(rhoVy)
      ,CHFp_FRA1(rhoVz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_REAL(vpx)
      ,CHFp_CONST_REAL(vpy)
      ,CHFp_CONST_REAL(vpz)
      ,CHFp_CONST_REAL(qp) );

#define FORT_TSC_DEPOSIT_CURRENT FORTRAN_NAME( inlineTSC_DEPOSIT_CURRENT, inlineTSC_DEPOSIT_CURRENT)
#define FORTNT_TSC_DEPOSIT_CURRENT FORTRAN_NAME( TSC_DEPOSIT_CURRENT, tsc_deposit_current)

inline void FORTRAN_NAME(inlineTSC_DEPOSIT_CURRENT, inlineTSC_DEPOSIT_CURRENT)(
      CHFp_FRA1(rhoVx)
      ,CHFp_FRA1(rhoVy)
      ,CHFp_FRA1(rhoVz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_REAL(vpx)
      ,CHFp_CONST_REAL(vpy)
      ,CHFp_CONST_REAL(vpz)
      ,CHFp_CONST_REAL(qp) )
{
 CH_TIMELEAF("FORT_TSC_DEPOSIT_CURRENT");
 FORTRAN_NAME( TSC_DEPOSIT_CURRENT ,tsc_deposit_current )(
      CHFt_FRA1(rhoVx)
      ,CHFt_FRA1(rhoVy)
      ,CHFt_FRA1(rhoVz)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp)
      ,CHFt_CONST_REAL(vpx)
      ,CHFt_CONST_REAL(vpy)
      ,CHFt_CONST_REAL(vpz)
      ,CHFt_CONST_REAL(qp) );
}
#endif  // GUARDTSC_DEPOSIT_CURRENT 

#ifndef GUARDCIC_INTERPOLATE_FIELDS 
#define GUARDCIC_INTERPOLATE_FIELDS 
// Prototype for Fortran procedure cic_interpolate_fields ...
//
void FORTRAN_NAME( CIC_INTERPOLATE_FIELDS ,cic_interpolate_fields )(
      CHFp_REAL(Epx)
      ,CHFp_REAL(Epy)
      ,CHFp_REAL(Epz)
      ,CHFp_REAL(Bpx)
      ,CHFp_REAL(Bpy)
      ,CHFp_REAL(Bpz)
      ,CHFp_CONST_FRA1(Ex)
      ,CHFp_CONST_FRA1(Ey)
      ,CHFp_CONST_FRA1(Ez)
      ,CHFp_CONST_FRA1(Bx)
      ,CHFp_CONST_FRA1(By)
      ,CHFp_CONST_FRA1(Bz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) );

#define FORT_CIC_INTERPOLATE_FIELDS FORTRAN_NAME( inlineCIC_INTERPOLATE_FIELDS, inlineCIC_INTERPOLATE_FIELDS)
#define FORTNT_CIC_INTERPOLATE_FIELDS FORTRAN_NAME( CIC_INTERPOLATE_FIELDS, cic_interpolate_fields)

inline void FORTRAN_NAME(inlineCIC_INTERPOLATE_FIELDS, inlineCIC_INTERPOLATE_FIELDS)(
      CHFp_REAL(Epx)
      ,CHFp_REAL(Epy)
      ,CHFp_REAL(Epz)
      ,CHFp_REAL(Bpx)
      ,CHFp_REAL(Bpy)
      ,CHFp_REAL(Bpz)
      ,CHFp_CONST_FRA1(Ex)
      ,CHFp_CONST_FRA1(Ey)
      ,CHFp_CONST_FRA1(Ez)
      ,CHFp_CONST_FRA1(Bx)
      ,CHFp_CONST_FRA1(By)
      ,CHFp_CONST_FRA1(Bz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) )
{
 CH_TIMELEAF("FORT_CIC_INTERPOLATE_FIELDS");
 FORTRAN_NAME( CIC_INTERPOLATE_FIELDS ,cic_interpolate_fields )(
      CHFt_REAL(Epx)
      ,CHFt_REAL(Epy)
      ,CHFt_REAL(Epz)
      ,CHFt_REAL(Bpx)
      ,CHFt_REAL(Bpy)
      ,CHFt_REAL(Bpz)
      ,CHFt_CONST_FRA1(Ex)
      ,CHFt_CONST_FRA1(Ey)
      ,CHFt_CONST_FRA1(Ez)
      ,CHFt_CONST_FRA1(Bx)
      ,CHFt_CONST_FRA1(By)
      ,CHFt_CONST_FRA1(Bz)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp) );
}
#endif  // GUARDCIC_INTERPOLATE_FIELDS 

#ifndef GUARDTSC_INTERPOLATE_FIELDS 
#define GUARDTSC_INTERPOLATE_FIELDS 
// Prototype for Fortran procedure tsc_interpolate_fields ...
//
void FORTRAN_NAME( TSC_INTERPOLATE_FIELDS ,tsc_interpolate_fields )(
      CHFp_REAL(Epx)
      ,CHFp_REAL(Epy)
      ,CHFp_REAL(Epz)
      ,CHFp_REAL(Bpx)
      ,CHFp_REAL(Bpy)
      ,CHFp_REAL(Bpz)
      ,CHFp_CONST_FRA1(Ex)
      ,CHFp_CONST_FRA1(Ey)
      ,CHFp_CONST_FRA1(Ez)
      ,CHFp_CONST_FRA1(Bx)
      ,CHFp_CONST_FRA1(By)
      ,CHFp_CONST_FRA1(Bz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) );

#define FORT_TSC_INTERPOLATE_FIELDS FORTRAN_NAME( inlineTSC_INTERPOLATE_FIELDS, inlineTSC_INTERPOLATE_FIELDS)
#define FORTNT_TSC_INTERPOLATE_FIELDS FORTRAN_NAME( TSC_INTERPOLATE_FIELDS, tsc_interpolate_fields)

inline void FORTRAN_NAME(inlineTSC_INTERPOLATE_FIELDS, inlineTSC_INTERPOLATE_FIELDS)(
      CHFp_REAL(Epx)
      ,CHFp_REAL(Epy)
      ,CHFp_REAL(Epz)
      ,CHFp_REAL(Bpx)
      ,CHFp_REAL(Bpy)
      ,CHFp_REAL(Bpz)
      ,CHFp_CONST_FRA1(Ex)
      ,CHFp_CONST_FRA1(Ey)
      ,CHFp_CONST_FRA1(Ez)
      ,CHFp_CONST_FRA1(Bx)
      ,CHFp_CONST_FRA1(By)
      ,CHFp_CONST_FRA1(Bz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) )
{
 CH_TIMELEAF("FORT_TSC_INTERPOLATE_FIELDS");
 FORTRAN_NAME( TSC_INTERPOLATE_FIELDS ,tsc_interpolate_fields )(
      CHFt_REAL(Epx)
      ,CHFt_REAL(Epy)
      ,CHFt_REAL(Epz)
      ,CHFt_REAL(Bpx)
      ,CHFt_REAL(Bpy)
      ,CHFt_REAL(Bpz)
      ,CHFt_CONST_FRA1(Ex)
      ,CHFt_CONST_FRA1(Ey)
      ,CHFt_CONST_FRA1(Ez)
      ,CHFt_CONST_FRA1(Bx)
      ,CHFt_CONST_FRA1(By)
      ,CHFt_CONST_FRA1(Bz)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp) );
}
#endif  // GUARDTSC_INTERPOLATE_FIELDS 

#ifndef GUARDCIC_INTERPOLATE_MAGNETIC_FIELD 
#define GUARDCIC_INTERPOLATE_MAGNETIC_FIELD 
// Prototype for Fortran procedure cic_interpolate_magnetic_field ...
//
void FORTRAN_NAME( CIC_INTERPOLATE_MAGNETIC_FIELD ,cic_interpolate_magnetic_field )(
      CHFp_REAL(Bpx)
      ,CHFp_REAL(Bpy)
      ,CHFp_REAL(Bpz)
      ,CHFp_CONST_FRA1(Bgx)
      ,CHFp_CONST_FRA1(Bgy)
      ,CHFp_CONST_FRA1(Bgz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) );

#define FORT_CIC_INTERPOLATE_MAGNETIC_FIELD FORTRAN_NAME( inlineCIC_INTERPOLATE_MAGNETIC_FIELD, inlineCIC_INTERPOLATE_MAGNETIC_FIELD)
#define FORTNT_CIC_INTERPOLATE_MAGNETIC_FIELD FORTRAN_NAME( CIC_INTERPOLATE_MAGNETIC_FIELD, cic_interpolate_magnetic_field)

inline void FORTRAN_NAME(inlineCIC_INTERPOLATE_MAGNETIC_FIELD, inlineCIC_INTERPOLATE_MAGNETIC_FIELD)(
      CHFp_REAL(Bpx)
      ,CHFp_REAL(Bpy)
      ,CHFp_REAL(Bpz)
      ,CHFp_CONST_FRA1(Bgx)
      ,CHFp_CONST_FRA1(Bgy)
      ,CHFp_CONST_FRA1(Bgz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) )
{
 CH_TIMELEAF("FORT_CIC_INTERPOLATE_MAGNETIC_FIELD");
 FORTRAN_NAME( CIC_INTERPOLATE_MAGNETIC_FIELD ,cic_interpolate_magnetic_field )(
      CHFt_REAL(Bpx)
      ,CHFt_REAL(Bpy)
      ,CHFt_REAL(Bpz)
      ,CHFt_CONST_FRA1(Bgx)
      ,CHFt_CONST_FRA1(Bgy)
      ,CHFt_CONST_FRA1(Bgz)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp) );
}
#endif  // GUARDCIC_INTERPOLATE_MAGNETIC_FIELD 

#ifndef GUARDCIC_INTERPOLATE_DIR 
#define GUARDCIC_INTERPOLATE_DIR 
// Prototype for Fortran procedure cic_interpolate_dir ...
//
void FORTRAN_NAME( CIC_INTERPOLATE_DIR ,cic_interpolate_dir )(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_INTVECT(stag) );

#define FORT_CIC_INTERPOLATE_DIR FORTRAN_NAME( inlineCIC_INTERPOLATE_DIR, inlineCIC_INTERPOLATE_DIR)
#define FORTNT_CIC_INTERPOLATE_DIR FORTRAN_NAME( CIC_INTERPOLATE_DIR, cic_interpolate_dir)

inline void FORTRAN_NAME(inlineCIC_INTERPOLATE_DIR, inlineCIC_INTERPOLATE_DIR)(
      CHFp_REAL(particle_field)
      ,CHFp_CONST_FRA1(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_INTVECT(stag) )
{
 CH_TIMELEAF("FORT_CIC_INTERPOLATE_DIR");
 FORTRAN_NAME( CIC_INTERPOLATE_DIR ,cic_interpolate_dir )(
      CHFt_REAL(particle_field)
      ,CHFt_CONST_FRA1(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp)
      ,CHFt_CONST_INTVECT(stag) );
}
#endif  // GUARDCIC_INTERPOLATE_DIR 

#ifndef GUARDTSC_DEPOSIT 
#define GUARDTSC_DEPOSIT 
// Prototype for Fortran procedure tsc_deposit ...
//
void FORTRAN_NAME( TSC_DEPOSIT ,tsc_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_INTVECT(stag)
      ,CHFp_CONST_REAL(kernal)
      ,CHFp_CONST_REAL(q) );

#define FORT_TSC_DEPOSIT FORTRAN_NAME( inlineTSC_DEPOSIT, inlineTSC_DEPOSIT)
#define FORTNT_TSC_DEPOSIT FORTRAN_NAME( TSC_DEPOSIT, tsc_deposit)

inline void FORTRAN_NAME(inlineTSC_DEPOSIT, inlineTSC_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp)
      ,CHFp_CONST_INTVECT(stag)
      ,CHFp_CONST_REAL(kernal)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_TSC_DEPOSIT");
 FORTRAN_NAME( TSC_DEPOSIT ,tsc_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp)
      ,CHFt_CONST_INTVECT(stag)
      ,CHFt_CONST_REAL(kernal)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDTSC_DEPOSIT 

#ifndef GUARDTSC_INTERPOLATE 
#define GUARDTSC_INTERPOLATE 
// Prototype for Fortran procedure tsc_interpolate ...
//
void FORTRAN_NAME( TSC_INTERPOLATE ,tsc_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_TSC_INTERPOLATE FORTRAN_NAME( inlineTSC_INTERPOLATE, inlineTSC_INTERPOLATE)
#define FORTNT_TSC_INTERPOLATE FORTRAN_NAME( TSC_INTERPOLATE, tsc_interpolate)

inline void FORTRAN_NAME(inlineTSC_INTERPOLATE, inlineTSC_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_TSC_INTERPOLATE");
 FORTRAN_NAME( TSC_INTERPOLATE ,tsc_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDTSC_INTERPOLATE 

#ifndef GUARDW4_DEPOSIT 
#define GUARDW4_DEPOSIT 
// Prototype for Fortran procedure w4_deposit ...
//
void FORTRAN_NAME( W4_DEPOSIT ,w4_deposit )(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) );

#define FORT_W4_DEPOSIT FORTRAN_NAME( inlineW4_DEPOSIT, inlineW4_DEPOSIT)
#define FORTNT_W4_DEPOSIT FORTRAN_NAME( W4_DEPOSIT, w4_deposit)

inline void FORTRAN_NAME(inlineW4_DEPOSIT, inlineW4_DEPOSIT)(
      CHFp_FRA1(rho)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x)
      ,CHFp_CONST_REAL(q) )
{
 CH_TIMELEAF("FORT_W4_DEPOSIT");
 FORTRAN_NAME( W4_DEPOSIT ,w4_deposit )(
      CHFt_FRA1(rho)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x)
      ,CHFt_CONST_REAL(q) );
}
#endif  // GUARDW4_DEPOSIT 

#ifndef GUARDW4_INTERPOLATE 
#define GUARDW4_INTERPOLATE 
// Prototype for Fortran procedure w4_interpolate ...
//
void FORTRAN_NAME( W4_INTERPOLATE ,w4_interpolate )(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) );

#define FORT_W4_INTERPOLATE FORTRAN_NAME( inlineW4_INTERPOLATE, inlineW4_INTERPOLATE)
#define FORTNT_W4_INTERPOLATE FORTRAN_NAME( W4_INTERPOLATE, w4_interpolate)

inline void FORTRAN_NAME(inlineW4_INTERPOLATE, inlineW4_INTERPOLATE)(
      CHFp_REALVECT(particle_field)
      ,CHFp_CONST_FRA(field)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(x) )
{
 CH_TIMELEAF("FORT_W4_INTERPOLATE");
 FORTRAN_NAME( W4_INTERPOLATE ,w4_interpolate )(
      CHFt_REALVECT(particle_field)
      ,CHFt_CONST_FRA(field)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(x) );
}
#endif  // GUARDW4_INTERPOLATE 

#ifndef GUARDCIC_FIELD_DERIVATIVES 
#define GUARDCIC_FIELD_DERIVATIVES 
// Prototype for Fortran procedure cic_field_derivatives ...
//
void FORTRAN_NAME( CIC_FIELD_DERIVATIVES ,cic_field_derivatives )(
      CHFp_REALVECT(dEpx)
      ,CHFp_REALVECT(dEpy)
      ,CHFp_REALVECT(dEpz)
      ,CHFp_REALVECT(dBpx)
      ,CHFp_REALVECT(dBpy)
      ,CHFp_REALVECT(dBpz)
      ,CHFp_CONST_FRA1(Ex)
      ,CHFp_CONST_FRA1(Ey)
      ,CHFp_CONST_FRA1(Ez)
      ,CHFp_CONST_FRA1(Bx)
      ,CHFp_CONST_FRA1(By)
      ,CHFp_CONST_FRA1(Bz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) );

#define FORT_CIC_FIELD_DERIVATIVES FORTRAN_NAME( inlineCIC_FIELD_DERIVATIVES, inlineCIC_FIELD_DERIVATIVES)
#define FORTNT_CIC_FIELD_DERIVATIVES FORTRAN_NAME( CIC_FIELD_DERIVATIVES, cic_field_derivatives)

inline void FORTRAN_NAME(inlineCIC_FIELD_DERIVATIVES, inlineCIC_FIELD_DERIVATIVES)(
      CHFp_REALVECT(dEpx)
      ,CHFp_REALVECT(dEpy)
      ,CHFp_REALVECT(dEpz)
      ,CHFp_REALVECT(dBpx)
      ,CHFp_REALVECT(dBpy)
      ,CHFp_REALVECT(dBpz)
      ,CHFp_CONST_FRA1(Ex)
      ,CHFp_CONST_FRA1(Ey)
      ,CHFp_CONST_FRA1(Ez)
      ,CHFp_CONST_FRA1(Bx)
      ,CHFp_CONST_FRA1(By)
      ,CHFp_CONST_FRA1(Bz)
      ,CHFp_CONST_REALVECT(left_edge)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(xp) )
{
 CH_TIMELEAF("FORT_CIC_FIELD_DERIVATIVES");
 FORTRAN_NAME( CIC_FIELD_DERIVATIVES ,cic_field_derivatives )(
      CHFt_REALVECT(dEpx)
      ,CHFt_REALVECT(dEpy)
      ,CHFt_REALVECT(dEpz)
      ,CHFt_REALVECT(dBpx)
      ,CHFt_REALVECT(dBpy)
      ,CHFt_REALVECT(dBpz)
      ,CHFt_CONST_FRA1(Ex)
      ,CHFt_CONST_FRA1(Ey)
      ,CHFt_CONST_FRA1(Ez)
      ,CHFt_CONST_FRA1(Bx)
      ,CHFt_CONST_FRA1(By)
      ,CHFt_CONST_FRA1(Bz)
      ,CHFt_CONST_REALVECT(left_edge)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(xp) );
}
#endif  // GUARDCIC_FIELD_DERIVATIVES 

}

#endif
