program main_attributes
!   $gfortran-5 -o error error.f90 main_error.f90 -I/usr/local/include -L/usr/lib/x86_64-linux-gnu -lfftw3_omp -lfftw3 -lm -fopenmp
!!  error.f90:37:0:

!!       DO indiceY=1+2*window,ny-2*window
!!   ^
!!  internal compiler error: in gfc_omp_clause_default_ctor, at fortran/trans-openmp.c:481
!!  Please submit a full bug report,
!!  with preprocessed source if appropriate.
!!  See <file:///usr/share/doc/gcc-5/README.Bugs> for instructions.
!!  *** gfortran  4.8.4 ***
!!gfortran -o error error.f90 main_error.f90 -I/usr/local/include -L/usr/lib/x86_64-linux-gnu -lfftw3_omp -lfftw3 -lm -fopenmp
!!error.f90: In function ‘error’:
!!error.f90:37:0: internal compiler error: in gfc_omp_clause_default_ctor, at fortran/trans-openmp.c:172
!!     DO indiceY=1+2*window,ny-2*window
!! ^
!!Please submit a full bug report,
!!with preprocessed source if appropriate.
!!See <file:///usr/share/doc/gcc-4.8/README.Bugs> for instructions.

USE ERROR_MAG
USE, intrinsic:: iso_c_binding
IMPLICIT NONE
INCLUDE 'fftw3.f'
INTEGER(kind=2),ALLOCATABLE  :: datajp(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE:: datajp_c(:,:,:)
REAL,ALLOCATABLE :: matrix_outX(:,:,:)
INTEGER::Nx=500,Ny=500,Nz=500,i

ALLOCATE(datajp(Nz,Ny,Nx),datajp_c(Nz,Ny,Nx))
ALLOCATE(matrix_outX(Nz,Ny,Nx))

datajp=reshape((/(i, i=1,Nx*Ny*Nz)/),shape(datajp))
datajp_c=DCMPLX(datajp)

CALL error(datajp_c, Nx, Ny, Nz, 5, matrix_outX)
write(*,*) 'MAIN-IMAGEOUT X sum main= ',SUM(matrix_outX)

OPEN(44, FILE="outputfile.raw", FORM="UNFORMATTED", STATUS="NEW", ACTION="WRITE", ACCESS='STREAM')
write(44) matrix_outX
CLOSE(44)

DEALLOCATE(matrix_outX, datajp, datajp_c)

end program main_attributes
