#include <define.h>

MODULE MOD_SpinupSpeedup

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!   Speed up spinup process.
!
!  Created by Shupeng Zhang, 04/2026
!----------------------------------------------------------------------------

   integer :: nyearadjust = 0

   PUBLIC :: spinup_speedup

CONTAINS

   !-----------------------------------------------------------------------
   SUBROUTINE spinup_speedup ( idate, deltim, is_spinup)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_TimeManager
   USE MOD_Grid
   USE MOD_Vars_Global,         only: WATERBODY, nl_lake, spval
   USE MOD_LandPatch,           only: numpatch,  landpatch
   USE MOD_Vars_TimeInvariants, only: patchtype, lakedepth, dz_lake
   USE MOD_Vars_TimeVariables,  only: wdsrf, t_lake, lake_icefrac
   USE MOD_CheckEquilibrium,    only: tws_last
   USE MOD_Forcing,             only: gforc
   USE MOD_HistGridded,         only: ghist
   USE MOD_VectorMapWrite,      only: map_and_write_vector
   USE MOD_Lake,                only: adjust_lake_layer

   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   logical,  intent(in) :: is_spinup

   ! Local variables
   logical :: do_adjust
   integer :: i
   type(grid_type)       :: gridcheck
   real(r8), allocatable :: vecadjust(:)
   logical,  allocatable :: filter   (:)
   character(len=256)    :: filename

#if (defined CatchLateralFlow)

      do_adjust = is_spinup .and. isendofyear (idate, deltim)

      IF (do_adjust) THEN

         nyearadjust = nyearadjust + 1

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               allocate (vecadjust (numpatch));   vecadjust = spval
               allocate (filter    (numpatch));   filter    = .true.

               DO i = 1, numpatch
                  IF (patchtype(i) == 4) THEN
                     vecadjust(i) = max(wdsrf(i)-lakedepth(i)*1.e3, 0.)
                  ELSE
                     vecadjust(i) = wdsrf(i)
                  ENDIF

                  wdsrf(i) = wdsrf(i) - vecadjust(i)

                  IF (patchtype(i) == 4) THEN
                     IF (vecadjust(i) > 0.) THEN
                        dz_lake(:,i) = dz_lake(:,i) * lakedepth(i)/sum(dz_lake(:,i))
                        CALL adjust_lake_layer (nl_lake, dz_lake(:,i), t_lake(:,i), lake_icefrac(:,i))
                     ENDIF
                  ENDIF

                  IF (DEF_CheckEquilibrium) THEN
                     tws_last(i) = tws_last(i) - vecadjust(i)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

         IF (DEF_HISTORY_IN_VECTOR) THEN
            CALL gridcheck%define_by_copy (gforc)
         ELSE
            CALL gridcheck%define_by_copy (ghist)
         ENDIF

         filename = trim(DEF_dir_history) // '/' // trim(DEF_CASE_NAME) // '_spinup.nc'

         CALL map_and_write_vector ( &
            gridcheck, landpatch, vecadjust, filter, filename, 'wdsrf_adjust', &
            'Adjustment in surface water depth', 'mm', nyearadjust, 'iyear')

         IF (allocated(vecadjust)) deallocate(vecadjust)
         IF (allocated(filter   )) deallocate(filter   )

      ENDIF

#endif

   END SUBROUTINE spinup_speedup

END MODULE MOD_SpinupSpeedup
