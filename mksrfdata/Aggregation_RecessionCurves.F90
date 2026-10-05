#include <define.h>

SUBROUTINE Aggregation_RecessionCurves ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Create recession curve for the modeling reolustion
!
!  Original author: Shupeng Zhang, 09/2026
!
! !REVISIONS:

!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_NetCDFBlock
   USE MOD_NetCDFVector
   USE MOD_AggregationRequestData
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_Utils
   USE MOD_UserDefFun
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   integer :: ipatch, L, np, ip
   integer :: wmo_src
   real(r8):: sumarea

   type (block_data_real8_2d) :: ths_grid        ! dimensionless, [0 - 1]
   type (block_data_real8_2d) :: c0_grid         ! [log(m/(s^2) * (m/s)^(-c1))]
   type (block_data_real8_2d) :: c1_grid         ! dimensionless
   type (block_data_real8_2d) :: qmax_grid       ! [m/s]

   real(r8), allocatable :: ths_rcc_patches  (:)  ! dimensionless, [0-1]
   real(r8), allocatable :: df_rcc_patches   (:)  ! [(m/s)^ep * 1/m]
   real(r8), allocatable :: ep_rcc_patches   (:)  ! dimensionless
   real(r8), allocatable :: qmax_rcc_patches (:)  ! [m/s]

   real(r8), allocatable :: ths_one  (:)         ! dimensionless, [0 - 1]
   real(r8), allocatable :: c0_one   (:)         ! [log(m/(s^2) * (m/s)^(-c1))]
   real(r8), allocatable :: c1_one   (:)         ! dimensionless
   real(r8), allocatable :: qmax_one (:)         ! [m/s]
   real(r8), allocatable :: area_one (:)
   logical,  allocatable :: filter   (:)

   ! local variables for estimating the upscaled parameters using the Levenberg-Marquardt fitting method
   integer, parameter :: npoint = 12
   real(r8),parameter :: xdat(npoint) = (/0., 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4/)
                                 ! points of soil depth used for fitting recession curves [unit: m]
   real(r8) :: ydat  (npoint)    ! mean of baseflow at fine grids
   real(r8) :: ythis (npoint)    ! baseflow at oen fine grid

   integer, parameter :: nv = 3  ! number of fitted parameters
   real(r8)           :: xv(nv)  ! parameters to be fitted

   integer :: info

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#endif


      ! ---------------------------------------------------------------------------------
      !  aggregate the parameters from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------------------------
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/recession_curve/' // trim(cyear)

      IF (p_is_master) THEN
         write(*,'(/, A41)') 'Aggregate Recesstion Curve Parameters ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate ( ths_rcc_patches  (numpatch) )
            allocate ( df_rcc_patches   (numpatch) )
            allocate ( ep_rcc_patches   (numpatch) )
            allocate ( qmax_rcc_patches (numpatch) )

            ths_rcc_patches  (:) = spval
            df_rcc_patches   (:) = spval
            ep_rcc_patches   (:) = spval
            qmax_rcc_patches (:) = spval
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         lndname = '/tera13/zhangsp/datasets/StreamFlow/recession_parameters_merged.nc'

         ! lndname = trim(dir_rawdata)//'/recession_curves.nc'
         CALL allocate_block_data (gland, ths_grid)
         CALL ncio_read_block (lndname, 'theta_s', gland, ths_grid)

         ! lndname = trim(dir_rawdata)//'/recession_curves.nc'
         CALL allocate_block_data (gland, c0_grid  )
         CALL ncio_read_block (lndname, 'c0',   gland, c0_grid  )

         ! lndname = trim(dir_rawdata)//'/recession_curves.nc'
         CALL allocate_block_data (gland, c1_grid  )
         CALL ncio_read_block (lndname, 'c1',   gland, c1_grid  )

         ! lndname = trim(dir_rawdata)//'/recession_curves.nc'
         CALL allocate_block_data (gland, qmax_grid)
         CALL ncio_read_block (lndname, 'qmax', gland, qmax_grid)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, &
            data_r8_2d_in1 = ths_grid, data_r8_2d_in2 = c0_grid, &
            data_r8_2d_in3 = c1_grid,  data_r8_2d_in4 = qmax_grid)
#endif
      ENDIF

      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
               wmo_src = wmo_source (landpatch%ielm(ipatch))

               ths_rcc_patches  (ipatch) = ths_rcc_patches  (wmo_src)
               df_rcc_patches   (ipatch) = df_rcc_patches   (wmo_src)
               ep_rcc_patches   (ipatch) = ep_rcc_patches   (wmo_src)
               qmax_rcc_patches (ipatch) = qmax_rcc_patches (wmo_src)

               CYCLE
            ENDIF

            L = landpatch%settyp(ipatch)

            IF (L /= 0) THEN

               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
                  area = area_one, &
                  data_r8_2d_in1 = ths_grid,  data_r8_2d_out1 = ths_one, &
                  data_r8_2d_in2 = c0_grid,   data_r8_2d_out2 = c0_one,  &
                  data_r8_2d_in3 = c1_grid,   data_r8_2d_out3 = c1_one,  &
                  data_r8_2d_in4 = qmax_grid, data_r8_2d_out4 = qmax_one )

               allocate (filter (size(area_one)))

               filter = .not. &
                  (      isnan_ud(c0_one)   .or. (c0_one   < -300.) &
                    .or. isnan_ud(c1_one)   .or. (c1_one   <= 0.  ) &
                    .or. isnan_ud(qmax_one) .or. (qmax_one <= 0.  ) &
                    .or. isnan_ud(ths_one)  .or. (ths_one  <= 0.  ) )

               np = count(filter)

               IF( np > 1 ) THEN

                  area_one(1:np) = pack(area_one, filter)
                  ths_one (1:np) = pack(ths_one , filter)
                  c0_one  (1:np) = pack(c0_one  , filter)
                  c1_one  (1:np) = pack(c1_one  , filter)
                  qmax_one(1:np) = pack(qmax_one, filter)

                  sumarea = sum(area_one(1:np))
                  ths_rcc_patches  (ipatch) = sum (ths_one(1:np) * (area_one(1:np)/sumarea))

                  ydat = 0
                  DO ip = 1,np
                     CALL recession_curve ( (/exp(c0_one(ip)), 2.-c1_one(ip), qmax_one(ip)/), &
                        ths_rcc_patches(ipatch), npoint, xdat, ythis)
                     ydat = ydat + ythis * area_one(ip)
                  ENDDO
                  ydat(:) = ydat(:) / sumarea

                  qmax_rcc_patches (ipatch) = sum (qmax_one(1:np) * (area_one(1:np)/sumarea))
                  qmax_rcc_patches (ipatch) = max(qmax_rcc_patches(ipatch), maxval(ydat))

                  ep_rcc_patches   (ipatch) = 0.
                  df_rcc_patches   (ipatch) = sum(xdat*(log(qmax_rcc_patches(ipatch))-log(max(ydat,1.e-50)))) &
                                              /sum(xdat**2) /ths_rcc_patches(ipatch)
                  df_rcc_patches   (ipatch) = max(df_rcc_patches(ipatch), 1.e-4)

                  ! Fitting the van Genuchten SW retention parameters
                  xv(1) = df_rcc_patches  (ipatch)
                  xv(2) = ep_rcc_patches  (ipatch)
                  xv(3) = qmax_rcc_patches(ipatch)

                  CALL lmder_recession_curve ( npoint, nv, xv, ths_rcc_patches(ipatch), xdat, ydat, info )

                  IF (info /= 0) THEN
                     IF ( xv(1) >= 1.e-4 .and. xv(2) > -10. .and. xv(2) <= 10. &
                        .and. xv(3) > 0.  .and. xv(3) <= 4.e-5  ) THEN
                        df_rcc_patches  (ipatch) = xv(1)
                        ep_rcc_patches  (ipatch) = xv(2)
                        qmax_rcc_patches(ipatch) = xv(3)
                     ENDIF
                  ENDIF

               ENDIF

               deallocate(area_one)
               deallocate(ths_one )
               deallocate(c0_one  )
               deallocate(c1_one  )
               deallocate(qmax_one)
               deallocate(filter  )

            ENDIF

         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('ths_rcc_patches   ', ths_rcc_patches,  spval)
      CALL check_vector_data ('df_rcc_patches    ', df_rcc_patches,   spval)
      CALL check_vector_data ('ep_rcc_patches    ', ep_rcc_patches,   spval)
      CALL check_vector_data ('qmax_rcc_patches  ', qmax_rcc_patches, spval)
#endif

      ! Write-out parameters
      lndname = trim(landdir)//'/ths_rcc_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'ths_rcc_patches', 'patch', &
           landpatch, ths_rcc_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/df_rcc_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'df_rcc_patches', 'patch', &
           landpatch, df_rcc_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/ep_rcc_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'ep_rcc_patches', 'patch', &
           landpatch, ep_rcc_patches, DEF_Srfdata_CompressLevel)

      lndname = trim(landdir)//'/qmax_rcc_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'qmax_rcc_patches', 'patch', &
           landpatch, qmax_rcc_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/recession_curve_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (ths_rcc_patches, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'porosity',    compress = 1, write_mode = 'one', create_mode=.true.)
      CALL srfdata_map_and_write (df_rcc_patches, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'decayfactor', compress = 1, write_mode = 'one')
      CALL srfdata_map_and_write (ep_rcc_patches, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'exponent',    compress = 1, write_mode = 'one')
      CALL srfdata_map_and_write (qmax_rcc_patches, landpatch%settyp, typpatch, m_patch2diag, &
         spval, lndname, 'qmax',        compress = 1, write_mode = 'one')
#endif

      ! Deallocate the allocatable array
      IF (p_is_worker) THEN

         IF (allocated(ths_rcc_patches ))   deallocate (ths_rcc_patches )
         IF (allocated(df_rcc_patches  ))   deallocate (df_rcc_patches  )
         IF (allocated(ep_rcc_patches  ))   deallocate (ep_rcc_patches  )
         IF (allocated(qmax_rcc_patches))   deallocate (qmax_rcc_patches)

      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif


CONTAINS

   ! ----------
   SUBROUTINE recession_curve (x, ths, npoint, xdat, ydat)

   IMPLICIT NONE

      integer  :: npoint
      real(r8) :: x(3), ths, xdat(npoint)
      real(r8) :: ydat(npoint)

      ! Local Variables
      real(r8) :: sc, ep, qmax
      real(r8) :: xtemp(npoint)

      sc   = x(1)
      ep   = x(2)
      qmax = x(3)

      IF (ep == 0) THEN
         ydat(:) = qmax * exp(-sc*ths*xdat(:))
      ELSE
         xtemp = max(-sc*ths*xdat(:)*ep+qmax**ep,1.e-14)
         ydat(:) = xtemp**(1/ep)
      ENDIF

   END SUBROUTINE

   ! ----------
   SUBROUTINE recession_curve_jacobian (x, ths, npoint, xdat, fjac)

   IMPLICIT NONE

      integer  :: npoint
      real(r8) :: x(3), ths, xdat(npoint)
      real(r8) :: fjac(npoint,3)

      ! Local Variables
      real(r8) :: sc, ep, qmax
      real(r8) :: xtemp(npoint)

      sc   = x(1)
      ep   = x(2)
      qmax = x(3)

      IF (ep == 0) THEN
         fjac(:,1) = qmax * exp(-sc*ths*xdat(:)) * (-ths*xdat(:))
         fjac(:,2) = qmax * exp(-sc*ths*xdat(:)) * 0.5*sc*ths*xdat(:)*(-sc*ths*xdat(:)+2*log(qmax))
         fjac(:,3) = qmax * exp(-sc*ths*xdat(:)) * 1./qmax
      ELSE
         xtemp = max(-sc*ths*xdat(:)*ep+qmax**ep,1.e-14)
         fjac(:,1) = xtemp**(1/ep) * (-ths*xdat(:))/xtemp
         fjac(:,2) = xtemp**(1/ep) * ((-sc*ths*xdat(:)+(qmax**ep)*log(qmax))/(ep*xtemp)-log(xtemp)/(ep*ep))
         fjac(:,3) = xtemp**(1/ep) * qmax**(ep-1)/xtemp
      ENDIF

   END SUBROUTINE

   !----------------------------------------------------
   SUBROUTINE lmder_recession_curve ( m, n, x, ths, xdat, ydat, info )

   !*******************************************************************************
   !
   !  LMDER minimizes M functions in N variables by the Levenberg-Marquardt method
   !  implemented for fitting the recession curve parameters.
   !
   !  Author:
   !    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
   !    FORTRAN90 version by John Burkardt.
   !    Modified by Nan Wei, 2019/01
   !    Modified by Shupeng Zhang, 2026/09
   !
   !  Reference:
   !
   !    Jorge More, Burton Garbow, Kenneth Hillstrom, User Guide for MINPACK-1,
   !    Technical Report ANL-80-74, Argonne National Laboratory, 1980.
   !
   !  Parameters:
   !
   !    real ( kind = 8 ) FJAC(M,N), an M by N array.  The upper
   !    N by N submatrix of FJAC contains an upper triangular matrix R with
   !    diagonal elements of nonincreasing magnitude such that
   !      P' * ( JAC' * JAC ) * P = R' * R,
   !    where P is a permutation matrix and JAC is the final calculated jacobian.
   !    Column J of P is column IPVT(J) of the identity matrix.  The lower
   !    trapezoidal part of FJAC contains information generated during
   !    the computation of R.
   !
   !    real ( kind = 8 ) FTOL.  Termination occurs when both the actual
   !    and predicted relative reductions in the sum of squares are at most FTOL.
   !    Therefore, FTOL measures the relative error desired in the sum of
   !    squares.  FTOL should be nonnegative.
   !
   !    real ( kind = 8 ) XTOL.  Termination occurs when the relative error
   !    between two consecutive iterates is at most XTOL.  XTOL should be
   !    nonnegative.
   !
   !    real ( kind = 8 ) GTOL.  Termination occurs when the cosine of the
   !    angle between FVEC and any column of the jacobian is at most GTOL in
   !    absolute value.  Therefore, GTOL measures the orthogonality desired
   !    between the function vector and the columns of the jacobian.  GTOL should
   !    be nonnegative.
   !
   !    integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
   !    calls to FCN is at least MAXFEV by the end of an iteration.
   !
   !    real ( kind = 8 ) FACTOR, determines the initial step bound.  This
   !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
   !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
   !    in the interval (0.1, 100) with 100 the recommended value.
   !
   !    integer ( kind = 4 ) INFO, error flag.
   !    INFO is set as follows:
   !    0, improper input parameters.
   !    1, both actual and predicted relative reductions in the sum of
   !       squares are at most FTOL.
   !    2, relative error between two consecutive iterates is at most XTOL.
   !    3, conditions for INFO = 1 and INFO = 2 both hold.
   !    4, the cosine of the angle between FVEC and any column of the jacobian
   !       is at most GTOL in absolute value.
   !    5, number of calls to FCN has reached MAXFEV.
   !    6, FTOL is too small.  No further reduction in the sum of squares
   !       is possible.
   !    7, XTOL is too small.  No further improvement in the approximate
   !       solution X is possible.
   !    8, GTOL is too small.  FVEC is orthogonal to the columns of the
   !       jacobian to machine precision.
   !
   !    integer ( kind = 4 ) IPVT(N), defines a permutation matrix P
   !    such that JAC*P = Q*R, where JAC is the final calculated jacobian, Q is
   !    orthogonal (not stored), and R is upper triangular with diagonal
   !    elements of nonincreasing magnitude.  Column J of P is column
   !    IPVT(J) of the identity matrix.
   !
   !    real ( kind = 8 ) QTF(N), contains the first N elements of Q'*FVEC.
   !
   IMPLICIT NONE

   ! Input
   integer ( kind = 4 ) m       ! the number of points
   integer ( kind = 4 ) n       ! the number of variables
   real ( kind = 8 ) ths        ! fraction of soil that is voids
   real ( kind = 8 ) xdat (m)   ! points of soil depth used for fitting recession curves
   real ( kind = 8 ) ydat (m)   ! mean of baseflow at fine grids

   ! Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
   !    estimate of the solution vector.  On output X contains the final
   !    estimate of the solution vector.
   real ( kind = 8 ) x(n)

   ! Local variables
   real ( kind = 8 ) actred
   real ( kind = 8 ) delta
   real ( kind = 8 ) diag(n)
   real ( kind = 8 ) dirder
   real ( kind = 8 ) epsmch
   real ( kind = 8 ) fjac(m,n)
   real ( kind = 8 ) fnorm
   real ( kind = 8 ) fnorm1
   real ( kind = 8 ) fvec(m)
   real ( kind = 8 ) gnorm
   integer ( kind = 4 ) info
   integer ( kind = 4 ) ipvt(n)
   integer ( kind = 4 ) iter
   integer ( kind = 4 ) j
   integer ( kind = 4 ) l
   integer ( kind = 4 ) nfev
   real ( kind = 8 ) par
   real ( kind = 8 ) pnorm
   real ( kind = 8 ) prered
   real ( kind = 8 ) qtf(n)
   real ( kind = 8 ) ratio
   real ( kind = 8 ) sum2
   real ( kind = 8 ) temp
   real ( kind = 8 ) temp1
   real ( kind = 8 ) temp2
   real ( kind = 8 ) wa1(n)
   real ( kind = 8 ) wa2(n)
   real ( kind = 8 ) wa3(n)
   real ( kind = 8 ) wa4(m)
   real ( kind = 8 ) xnorm

   real(r8),parameter :: factor = 0.1
   real(r8),parameter :: ftol   = 1.0e-5
   real(r8),parameter :: xtol   = 1.0e-4
   real(r8),parameter :: gtol   = 0.0
   integer, parameter :: maxfev = 400
   logical, parameter :: pivot  = .true.


      epsmch = epsilon ( epsmch )

      info = 0
      nfev = 0
      !
      !  Check the input parameters for errors.
      !
      IF ( n <= 0 ) THEN
         go to 300
      ENDIF

      IF ( m < n ) THEN
         go to 300
      ENDIF

      !
      !  Evaluate the residual at the starting point and calculate its norm.
      !
      CALL recession_curve (x, ths, m, xdat, fvec)
      fvec = fvec - ydat
      nfev = nfev + 1

      fnorm = sqrt(sum(fvec(:)**2))
      !
      !  Initialize Levenberg-Marquardt parameter and iteration counter.
      !
      par = 0.0D+00
      iter = 1
      !
      !  Beginning of the outer loop.
      !
      DO
         !
         !  Calculate the jacobian matrix.
         !
         CALL recession_curve_jacobian (x, ths, m, xdat, fjac)

         !
         !     Compute the QR factorization of the jacobian.
         !
         CALL qrfac ( m, n, fjac, m, pivot, ipvt, n, wa1, wa2 )

         !     On the first iteration, scale according
         !     to the norms of the columns of the initial jacobian.
         !
         IF ( iter == 1 ) THEN

            diag(1:n) = wa2(1:n)
            DO j = 1, n
               IF ( wa2(j) == 0.0D+00 ) THEN
                  diag(j) = 1.0D+00
               ENDIF
            ENDDO
            !
            !        On the first iteration, calculate the norm of the scaled X
            !        and initialize the step bound DELTA.
            !
            wa3(1:n) = diag(1:n) * x(1:n)

            xnorm = enorm ( n, wa3 )

            IF ( xnorm == 0.0D+00 ) THEN
               delta = factor
            ELSE
               delta = factor * xnorm
            ENDIF

         ENDIF
         !
         !     Form Q'*FVEC and store the first N components in QTF.
         !
         wa4(1:m) = fvec(1:m)

         DO j = 1, n

            IF ( fjac(j,j) /= 0.0D+00 ) THEN
               sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
               temp = - sum2 / fjac(j,j)
               wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
            ENDIF

            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)

         ENDDO
         !
         !     Compute the norm of the scaled gradient.
         !
         gnorm = 0.0D+00

         IF ( fnorm /= 0.0D+00 ) THEN

            DO j = 1, n
               l = ipvt(j)
               IF ( wa2(l) /= 0.0D+00 ) THEN
                  sum2 = dot_product ( qtf(1:j), fjac(1:j,j) ) / fnorm
                  gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
               ENDIF
            ENDDO

         ENDIF
         !
         !     Test for convergence of the gradient norm.
         !
         IF ( gnorm <= gtol ) THEN
            info = 4
            go to 300
         ENDIF
         !
         !     Rescale if necessary.
         !
         DO j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
         ENDDO
         !
         !     Beginning of the inner loop.
         !
         DO
            !
            !     Determine the Levenberg-Marquardt parameter.

            CALL lmpar ( n, fjac, m, ipvt, diag, qtf, delta, par, wa1, wa2 )

            !        Store the direction p and x + p. calculate the norm of p.
            !
            wa1(1:n) = - wa1(1:n)
            wa2(1:n) = x(1:n) + wa1(1:n)
            wa3(1:n) = diag(1:n) * wa1(1:n)

            pnorm = enorm ( n, wa3 )
            !
            !        On the first iteration, adjust the initial step bound.
            !
            IF ( iter == 1 ) THEN
               delta = min ( delta, pnorm )
            ENDIF
            !
            !        Evaluate the function at x + p and calculate its norm.
            !
            CALL recession_curve (wa2, ths, m, xdat, wa4)
            wa4 = wa4 - ydat
            nfev = nfev + 1

            fnorm1 = enorm ( m, wa4 )
            !
            !        Compute the scaled actual reduction.
            !
            IF ( 0.1D+00 * fnorm1 < fnorm ) THEN
               actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
            ELSE
               actred = - 1.0D+00
            ENDIF
            !
            !        Compute the scaled predicted reduction and
            !        the scaled directional derivative.
            !
            DO j = 1, n
               wa3(j) = 0.0D+00
               l = ipvt(j)
               temp = wa1(l)
               wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
            ENDDO

            temp1 = enorm ( n, wa3 ) / fnorm
            temp2 = ( sqrt ( par ) * pnorm ) / fnorm
            prered = temp1 ** 2 + temp2 ** 2 / 0.5D+00
            dirder = - ( temp1 ** 2 + temp2 ** 2 )
            !
            !        Compute the ratio of the actual to the predicted reduction.
            !
            IF ( prered /= 0.0D+00 ) THEN
               ratio = actred / prered
            ELSE
               ratio = 0.0D+00
            ENDIF
            !
            !        Update the step bound.
            !
            IF ( ratio <= 0.25D+00 ) THEN

               IF ( 0.0D+00 <= actred ) THEN
                  temp = 0.5D+00
               ENDIF

               IF ( actred < 0.0D+00 ) THEN
                  temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
               ENDIF

               IF ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) THEN
                  temp = 0.1D+00
               ENDIF

               delta = temp * min ( delta, pnorm / 0.1D+00 )
               par = par / temp

            ELSE

               IF ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) THEN
                  delta = 2.0D+00 * pnorm
                  par = 0.5D+00 * par
               ENDIF

            ENDIF
            !
            !        Successful iteration.
            !
            !        Update X, FVEC, and their norms.
            !
            IF ( 0.0001D+00 <= ratio ) THEN
               x(1:n) = wa2(1:n)
               wa2(1:n) = diag(1:n) * x(1:n)
               fvec(1:m) = wa4(1:m)
               xnorm = enorm ( n, wa2 )
               fnorm = fnorm1
               iter = iter + 1
            ENDIF
            !
            !        Tests for convergence.
            !
            IF ( abs ( actred) <= ftol .and. &
               prered <= ftol .and. &
               0.5D+00 * ratio <= 1.0D+00 ) THEN
               info = 1
            ENDIF

            IF ( delta <= xtol * xnorm ) THEN
               info = 2
            ENDIF

            IF ( abs ( actred) <= ftol .and. prered <= ftol &
               .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) THEN
               info = 3
            ENDIF

            IF ( info /= 0 ) THEN
               go to 300
            ENDIF
            !
            !        Tests for termination and stringent tolerances.
            !
            IF ( nfev >= maxfev ) THEN
               info = 5
            ENDIF

            IF ( abs ( actred ) <= epsmch .and. prered <= epsmch &
               .and. 0.5D+00 * ratio <= 1.0D+00 ) THEN
               info = 6
            ENDIF

            IF ( delta <= epsmch * xnorm ) THEN
               info = 7
            ENDIF

            IF ( gnorm <= epsmch ) THEN
               info = 8
            ENDIF

            IF ( info /= 0 ) THEN
               go to 300
            ENDIF
            !
            !        End of the inner loop. repeat IF iteration unsuccessful.
            !
            IF ( 0.0001D+00 <= ratio ) THEN
               EXIT
            ENDIF

         ENDDO
         !
         !  End of the outer loop.
         !
      ENDDO

      300 continue

   END SUBROUTINE lmder_recession_curve

END SUBROUTINE Aggregation_RecessionCurves
!EOP
