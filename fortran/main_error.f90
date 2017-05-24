program main_attributes
!  gfortran -o main main_attributes.f90 dissimilarity_error_magnitude.f90 -I/usr/local/include -L/usr/lib/x86_64-linux-gnu 
! -lfftw3_omp -lfftw3 -lm -fopenmp
USE ERROR_MAG
USE, intrinsic:: iso_c_binding
IMPLICIT NONE
INCLUDE 'fftw3.f'
INTEGER(kind=2),ALLOCATABLE  :: datajp(:,:,:)
DOUBLE COMPLEX, ALLOCATABLE:: datajp_c(:,:,:)
REAL,ALLOCATABLE :: matrix_outX(:,:,:)
INTEGER::Nx,Ny,Nz,i

Nx = 500
Ny = 500
Nz = 500

ALLOCATE(datajp(Nz,Ny,Nx),datajp_c(Nz,Ny,Nx))
ALLOCATE(matrix_outX(Nz,Ny,Nx))

datajp=reshape((/(i, i=1,Nx*Ny*Nz)/),shape(datajp))
datajp_c=DCMPLX(datajp)

CALL error(datajp_c, 5, matrix_outX)
write(*,*) 'MAIN-IMAGEOUT X sum main= ',SUM(matrix_outX)

OPEN(44, FILE="outputfile.raw", FORM="UNFORMATTED", STATUS="NEW", ACTION="WRITE", ACCESS='STREAM')
write(44) matrix_outX
CLOSE(44)

DEALLOCATE(matrix_outX, datajp, datajp_c)

end program main_attributes
