#include <define.h>

MODULE MOD_VectorMapWrite

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!   Map vector to gridded data and write out.
!
!  Created by Shupeng Zhang, 04/2026
!----------------------------------------------------------------------------

CONTAINS

   SUBROUTINE map_and_write_vector ( &
         grid, pixelset, vector, filter, filename, varname, longname, units, &
         itime_in_file, timedimname, amount_in_grid)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_Block
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_NetCDFSerial
   USE MOD_SpatialMapping
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   type(grid_type),     intent(in) :: grid
   type(pixelset_type), intent(in) :: pixelset

   real(r8), intent(in) :: vector(:)
   logical,  intent(in) :: filter(:)

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: varname
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: timedimname

   logical, intent(in), optional :: amount_in_grid

   ! Local variables
   type(grid_concat_type)     :: gridconcat
   type(block_data_real8_2d)  :: sumarea, data_xy_2d
   type(spatial_mapping_type) :: mapping
   integer :: mpi_data_id = 3088
   integer :: xblk, yblk, xloc, yloc, xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
   integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
   integer :: rmesg(3), smesg(3), isrc
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)
   logical :: amount
   logical :: fexists


      CALL mapping%build_arealweighted (grid, pixelset)
      CALL gridconcat%set (grid)

      IF (p_is_io) CALL allocate_block_data (grid, sumarea)
      CALL mapping%get_sumarea (sumarea, filter)

      IF (p_is_io) CALL allocate_block_data (grid, data_xy_2d)
      CALL mapping%pset2grid (vector, data_xy_2d, spv = spval, msk = filter)

      amount = .false.
      IF (present(amount_in_grid)) amount = amount_in_grid

      IF (.not. amount) THEN
         IF (p_is_io) THEN
            DO iblkme = 1, gblock%nblkme
               xblk = gblock%xblkme(iblkme)
               yblk = gblock%yblkme(iblkme)

               DO yloc = 1, grid%ycnt(yblk)
                  DO xloc = 1, grid%xcnt(xblk)

                     IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                        IF (data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                           data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                              = data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                              / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                        ENDIF
                     ELSE
                        data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                     ENDIF

                  ENDDO
               ENDDO

            ENDDO
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         allocate (vdata (gridconcat%ginfo%nlon, gridconcat%ginfo%nlat))
         vdata(:,:) = spval
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN
         DO idata = 1, gridconcat%ndatablk

            CALL mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_data_id, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            ixseg = rmesg(2)
            iyseg = rmesg(3)

            xgdsp = gridconcat%xsegs(ixseg)%gdsp
            ygdsp = gridconcat%ysegs(iyseg)%gdsp
            xcnt  = gridconcat%xsegs(ixseg)%cnt
            ycnt  = gridconcat%ysegs(iyseg)%cnt

            allocate (rbuf(xcnt,ycnt))

            CALL mpi_recv (rbuf, xcnt*ycnt, MPI_REAL8, &
               isrc, mpi_data_id, p_comm_glb, p_stat, p_err)

            vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = rbuf
            deallocate (rbuf)

         ENDDO
      ENDIF

      IF (p_is_io) THEN
         DO iyseg = 1, gridconcat%nyseg
            DO ixseg = 1, gridconcat%nxseg

               iblk = gridconcat%xsegs(ixseg)%blk
               jblk = gridconcat%ysegs(iyseg)%blk

               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                  xbdsp = gridconcat%xsegs(ixseg)%bdsp
                  ybdsp = gridconcat%ysegs(iyseg)%bdsp
                  xcnt  = gridconcat%xsegs(ixseg)%cnt
                  ycnt  = gridconcat%ysegs(iyseg)%cnt

                  allocate (sbuf (xcnt,ycnt))
                  sbuf = data_xy_2d%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                  smesg = (/p_iam_glb, ixseg, iyseg/)
                  CALL mpi_send (smesg, 3, MPI_INTEGER, &
                     p_address_master, mpi_data_id, p_comm_glb, p_err)
                  CALL mpi_send (sbuf, xcnt*ycnt, MPI_REAL8, &
                     p_address_master, mpi_data_id, p_comm_glb, p_err)

                  deallocate (sbuf)

               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

#else
      DO iyseg = 1, gridconcat%nyseg
         DO ixseg = 1, gridconcat%nxseg
            iblk = gridconcat%xsegs(ixseg)%blk
            jblk = gridconcat%ysegs(iyseg)%blk
            IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
               xbdsp = gridconcat%xsegs(ixseg)%bdsp
               ybdsp = gridconcat%ysegs(iyseg)%bdsp
               xgdsp = gridconcat%xsegs(ixseg)%gdsp
               ygdsp = gridconcat%ysegs(iyseg)%gdsp
               xcnt  = gridconcat%xsegs(ixseg)%cnt
               ycnt  = gridconcat%ysegs(iyseg)%cnt

               vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                  data_xy_2d%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
            ENDIF
         ENDDO
      ENDDO
#endif

      IF (p_is_master) THEN

         inquire (file=trim(filename), exist=fexists)
         IF (.not. fexists) THEN
            CALL ncio_create_file (trim(filename))
         ENDIF

         IF (itime_in_file > 0) THEN
            IF (.not. ncio_var_exist (filename, timedimname, readflag = .false.)) THEN
               CALL ncio_define_dimension(filename, timedimname, 0)
            ENDIF
         ENDIF

         IF (.not. ncio_var_exist (filename, 'lat', readflag = .false.)) THEN
            CALL ncio_define_dimension(filename, 'lat', gridconcat%ginfo%nlat)
            CALL ncio_write_serial    (filename, 'lat', gridconcat%ginfo%lat_c, 'lat')
            CALL ncio_put_attr        (filename, 'lat', 'long_name', 'latitude')
            CALL ncio_put_attr        (filename, 'lat', 'units', 'degrees_north')
         ENDIF

         IF (.not. ncio_var_exist (filename, 'lon', readflag = .false.)) THEN
            CALL ncio_define_dimension(filename, 'lon', gridconcat%ginfo%nlon)
            CALL ncio_write_serial    (filename, 'lon', gridconcat%ginfo%lon_c, 'lon')
            CALL ncio_put_attr        (filename, 'lon', 'long_name', 'longitude')
            CALL ncio_put_attr        (filename, 'lon', 'units', 'degrees_east')
         ENDIF

         IF (itime_in_file >= 1) THEN
            CALL ncio_write_serial_time (filename, varname, itime_in_file, vdata, &
               'lon', 'lat', trim(timedimname), compress = 1)
         ELSE
            CALL ncio_write_serial (filename, varname, vdata, 'lon', 'lat', compress = 1)
         ENDIF

         IF (itime_in_file <= 1) THEN
            CALL ncio_put_attr (filename, varname, 'long_name', longname)
            CALL ncio_put_attr (filename, varname, 'units', units)
            CALL ncio_put_attr (filename, varname, 'missing_value', spval)
         ENDIF

         deallocate (vdata)

      ENDIF

   END SUBROUTINE map_and_write_vector

END MODULE MOD_VectorMapWrite
