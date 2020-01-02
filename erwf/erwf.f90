!==============================================================================
!>  Program to create hydrogenic radial wave-functions.
!
PROGRAM erwf

!    Written by A. Senchuk and C. Froese-Fischer        December 2019
!
!=============================================================================


  !> Print startup message
  WRITE (6, *) 'ERWF'
  WRITE (6, *) 'This program estimates radial wave functions'
  WRITE (6, *) 'for orbitals'
  WRITE (6, *) 'Input files: isodata, rcsf.inp, optional rwfn file'
  WRITE (6, *) 'Output file: rwfn.inp'



  !> Open rcsf.inp, isodata, and rwfn.inp (if present)

  
  !> Load the orbitals and nuclear data
  CALL load_orbs
  CALL load_isodata

  
  !> Generates radial grid 
  call radial_grid
  
  
  !> Generate radial wave functions
  CALL genrwf

  
  !> Orthogonalize the wave functions

  
  !> Write radial wavefunctions to file

  
  !> Perform all final tasks needed to end rwfnestimate
  

END PROGRAM erwf
