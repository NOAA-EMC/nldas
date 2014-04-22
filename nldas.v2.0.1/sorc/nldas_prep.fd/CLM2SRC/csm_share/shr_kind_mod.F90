!===============================================================================
! CVS: $Id: shr_kind_mod.F90,v 1.1.1.1 2003/02/06 16:10:55 jgottsch Exp $
! CVS: $Source: /YUKON/LDASCVS/LDAS/SRC/CLM2SRC/csm_share/shr_kind_mod.F90,v $
! CVS: $Name:  $
!===============================================================================

MODULE shr_kind_mod

   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
   integer,parameter :: SHR_KIND_R8 = selected_real_kind( 5) ! 8 byte real
   integer,parameter :: SHR_KIND_R4 = selected_real_kind( 5) ! 4 byte real
   integer,parameter :: SHR_KIND_RN = kind(1.0)              ! native real
   integer,parameter :: SHR_KIND_I8 = selected_int_kind ( 5) ! 8 byte integer
   integer,parameter :: SHR_KIND_I4 = selected_int_kind ( 5) ! 4 byte integer
   integer,parameter :: SHR_KIND_IN = kind(1)                ! native integer

END MODULE shr_kind_mod
