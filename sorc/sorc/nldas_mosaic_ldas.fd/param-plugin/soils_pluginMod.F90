!-------------------------------------------------------------------------
! NASA Goddard Space Flight Center Land Information System (LIS) V4.0.2
! Released October 2005
!
! See SOFTWARE DISTRIBUTION POLICY for software distribution policies
!
! The LIS source code and documentation are in the public domain,
! available without fee for educational, research, non-commercial and
! commercial purposes.  Users may distribute the binary or source
! code to third parties provided this statement appears on all copies and
! that no charge is made for such copies.
!
! NASA GSFC MAKES NO REPRESENTATIONS ABOUT THE SUITABILITY OF THE
! SOFTWARE FOR ANY PURPOSE.  IT IS PROVIDED AS IS WITHOUT EXPRESS OR
! IMPLIED WARRANTY.  NEITHER NASA GSFC NOR THE US GOVERNMENT SHALL BE
! LIABLE FOR ANY DAMAGES SUFFERED BY THE USER OF THIS SOFTWARE.
!
! See COPYRIGHT.TXT for copyright details.
!
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: soils_pluginMod.F90  
! 
! !DESCRIPTION: 
!   This module contains the definition of the functions used for
!   incorporating a new model forcing scheme. 
!   
! !REVISION HISTORY: 
!  11 Dec 03    Sujay Kumar  Initial Specification
! 
! !INTERFACE:
module soils_pluginMod
!EOP  
  implicit none
  
contains
!BOP
! !ROUTINE: soils_plugin
!
! !DESCRIPTION:
!
! This is a custom-defined plugin point for introducing a new forcing scheme.
! The interface mandates that the following routines be implemented
! and registered for each model forcing scheme. 
!
!  \begin{description}
!  \item[retrieval of forcing data]      
!      Routines to retrieve forcing data and to interpolate them. 
!      (to be registered using registerget)
!  \item[definition of native domain]
!      Routines to define the native domain as a kgds array
!      (to be registered using registerdefnat)
!  \item[temporal interpolation] 
!      Interpolate forcing data temporally. 
!      (to be registered using registertimeinterp)
!  \end{description}
! Multiple forcing schemes can be 
! included as well, each distinguished in the function table registry
! by the associated forcing index assigned in the card file. 
! 
! !INTERFACE:
  subroutine soils_plugin
    external read_faosand, read_faoclay, read_faosilt
    external read_statsgosand, read_statsgoclay, read_statsgosilt
    external read_gswpsand, read_gswpclay, read_gswpsilt
    external read_gswp_w_sat, read_gswp_w_sat_matp,    &
             read_gswp_w_sat_hydc, read_gswp_w_bpower, &
             read_gswp_w_wilt
    external read_gswp_soilclass
    external read_nldas_soilclass
! !USES:
    ! lis%d%soil
    call registerreadsand(2,read_faosand)    
    call registerreadclay(2,read_faoclay)
    call registerreadsilt(2,read_faosilt)
    
    call registerreadsand(3,read_statsgosand)    
    call registerreadclay(3,read_statsgoclay)
    call registerreadsilt(3,read_statsgosilt)

    call registerreadsand(4,read_gswpsand)    
    call registerreadclay(4,read_gswpclay)
    call registerreadsilt(4,read_gswpsilt)

    call registerreadsoilclass(5, read_gswp_soilclass)
    
    call registerreadsoilclass(6, read_nldas_soilclass)


    ! lis%p%soilp_type
    call registerreadwsat(2, read_gswp_w_sat)
    call registerreadwsatmatp(2, read_gswp_w_sat_matp)
    call registerreadwsathydc(2, read_gswp_w_sat_hydc)
    call registerreadwbpower(2, read_gswp_w_bpower)
    call registerreadwwilt(2, read_gswp_w_wilt)
!EOC
  end subroutine soils_plugin
end module soils_pluginMod
