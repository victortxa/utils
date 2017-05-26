MODULE ERROR_MAG
USE omp_lib
USE, intrinsic:: iso_c_binding
CONTAINS

SUBROUTINE error(matrix_in, nx, ny, nz, filter_window, matrix_outX)
IMPLICIT NONE
INCLUDE 'fftw3.f'
INTEGER, INTENT(IN):: filter_window, nx, ny, nz
DOUBLE COMPLEX, INTENT(IN):: matrix_in(nx,ny,nz)
REAL, INTENT(OUT), ALLOCATABLE:: matrix_outX(:, :, :)
DOUBLE COMPLEX, ALLOCATABLE:: OUTPUT_FFTW(:, :, :)
INTEGER:: indiceX, indiceY, indiceZ, window, iret, nthreads
integer*8:: plan

! FFTW output memory allocation
ALLOCATE(OUTPUT_FFTW(filter_window, filter_window, filter_window))

ALLOCATE(matrix_outX(nx, ny, nz))
matrix_outX = 0.
window = (filter_window-1)/2    
!***********************************
call dfftw_init_threads(iret)
call dfftw_plan_with_nthreads(nthreads)
!***********************************
CALL dfftw_plan_dft_3d(plan, filter_window, filter_window, filter_window, OUTPUT_FFTW, OUTPUT_FFTW, FFTW_FORWARD, FFTW_ESTIMATE)
!$OMP PARALLEL DO DEFAULT(SHARED) SHARED(matrix_in, window, filter_window, plan, nx, ny, nz) &
!$OMP PRIVATE(indiceX, indiceY, indiceZ, OUTPUT_FFTW, matrix_outX)
DO indiceZ=1+2*window,nz-2*window
    DO indiceY=1+2*window,ny-2*window
        DO indiceX=1+2*window,nx-2*window
            OUTPUT_FFTW = ABS(matrix_in(indiceX-2*window:indiceX, indiceY-window:indiceY+window,& 
                                                                 indiceZ-window:indiceZ+window) - &
                              matrix_in(indiceX:indiceX+2*window, indiceY-window:indiceY+window,&
                                                                 indiceZ-window:indiceZ+window) )
            CALL dfftw_execute_dft(plan, OUTPUT_FFTW, OUTPUT_FFTW)
            matrix_outX(indiceX, indiceY, indiceZ) = SUM(ABS(OUTPUT_FFTW))/filter_window**3
        END DO
    END DO
END DO
!$OMP END PARALLEL DO

CALL dfftw_destroy_plan(plan)
!***********************************
CALL dfftw_cleanup_threads()

DEALLOCATE(OUTPUT_FFTW)

END SUBROUTINE error
END MODULE ERROR_MAG
