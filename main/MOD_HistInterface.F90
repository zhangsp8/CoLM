#include <define.h>

MODULE MOD_HistInterface

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!     Write out gridded model results to history files.
!
!  Original version: Yongjiu Dai, September 15, 1999, 03/2014
!
! !REVISIONS:
!  Shupeng Zhang, 05/2023: 1) porting codes to MPI parallel version
!
!  TODO...(need complement)
!----------------------------------------------------------------------------

   USE MOD_Vars_1DAccFluxes
   USE MOD_Vars_Global, only: spval
   USE MOD_NetCDFSerial

   USE MOD_HistGridded
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE MOD_HistVector
#endif
#ifdef SinglePoint
   USE MOD_HistSingle
#endif

   USE MOD_Vars_1DAccFluxes, only: nac

   character(len=10) :: HistForm ! 'Gridded', 'Vector', 'Single'

   PUBLIC :: write_history_variable_2d
   PUBLIC :: write_history_variable_urb_2d
   PUBLIC :: write_history_variable_3d
   PUBLIC :: write_history_variable_4d
   PUBLIC :: hist_write_time

!--------------------------------------------------------------------------
CONTAINS

   SUBROUTINE write_history_variable_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units, acc_num, is_inst, average_mode)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8),         intent(inout) :: acc_vec(:)
   character(len=*), intent(in)    :: file_hist
   character(len=*), intent(in)    :: varname
   integer,          intent(in)    :: itime_in_file
   character(len=*), intent(in)    :: longname
   character(len=*), intent(in)    :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   real(r8), intent(in), optional  :: acc_num(:)
   logical,  intent(in), optional  :: is_inst

   character(len=*), intent(in), optional :: average_mode

   logical :: is_inst_

      IF (.not. is_hist) RETURN

      is_inst_ = .false.
      IF (present(is_inst)) is_inst_ = is_inst

      IF (.not. is_inst_) THEN
         IF ( .not. present(acc_num) ) THEN
            IF (p_is_worker) &
               WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
         ELSE
            IF (p_is_worker)  &
               WHERE (acc_vec/=spval .and. acc_num>0) acc_vec = acc_vec / acc_num
         ENDIF
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         IF (present(average_mode)) THEN
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units, average_mode)
         ELSE
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
         ENDIF
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         IF (present(average_mode)) THEN
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units, average_mode)
         ELSE
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
         ENDIF
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_2d


#ifdef URBAN_MODEL
   SUBROUTINE write_history_variable_urb_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units, acc_num)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   real(r8), intent(in), optional  :: acc_num(:)

      IF (.not. is_hist) RETURN

      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker)  &
            WHERE (acc_vec/=spval .and. acc_num>0) acc_vec = acc_vec / acc_num
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_urb_2d ( &
            acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         !TODO: currently, it is not applicable to urban variables
         CALL aggregate_to_vector_and_write_urb_2d ( &
            acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_urb_2d
#endif


   SUBROUTINE write_history_variable_3d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units

   ! Local variables
   integer :: iblkme, xblk, yblk, xloc, yloc, i1
   integer :: compress

      IF (.not. is_hist) RETURN

      IF (p_is_worker)  THEN
         WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_3d (acc_vec, file_hist, varname, itime_in_file, &
            dim1name, ndim1, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_3d


   SUBROUTINE write_history_variable_4d ( is_hist,   &
         acc_vec, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, dim2name, lb2, ndim2, &
         sumarea, filter, longname, units, acc_num)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:,:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name, dim2name
   integer, intent(in) :: lb1, ndim1, lb2, ndim2

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units
   real(r8), intent(in), optional:: acc_num(:)

   ! Local variables
   integer :: i1, i2

      IF (.not. is_hist) RETURN

      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker) THEN
            DO i1 = lbound(acc_vec,1), ubound(acc_vec,1)
               DO i2 = lbound(acc_vec,2), ubound(acc_vec,2)
                  WHERE (acc_vec(i1,i2,:)/=spval .and. acc_num>0) &
                        acc_vec(i1,i2,:) = acc_vec(i1,i2,:) / acc_num
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_4d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            dim2name, lb2, ndim2, sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_4d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            dim2name, lb2, ndim2, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_4d (acc_vec, file_hist, varname, itime_in_file, &
            dim1name, ndim1, dim2name, ndim2, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_4d


   SUBROUTINE hist_write_time (filename, filelast, dataname, time, itime)

   IMPLICIT NONE

   character (len=*), intent(in) :: filename
   character (len=*), intent(in) :: filelast
   character (len=*), intent(in) :: dataname
   integer, intent(in)  :: time(3)
   integer, intent(out) :: itime

      select CASE (HistForm)
      CASE ('Gridded')
         CALL hist_gridded_write_time (filename, filelast, dataname, time, itime)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL hist_vector_write_time  (filename, filelast, dataname, time, itime)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL hist_single_write_time  (filename, filelast, dataname, time, itime)
#endif
      END select

   END SUBROUTINE hist_write_time

END MODULE MOD_HistInterface
