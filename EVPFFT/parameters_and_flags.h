#pragma once

// make sure only on matrix_inverse macro is defined
#undef  LU_MATRIX_INVERSE
#define GJE_MATRIX_INVERSE

#define NON_SCHMID_EFFECTS
#define TWO_SIGN_SLIP_SYSTEMS

#define NPHMX  2       ! MAXIMUM # OF PHASES
#define NMODMX 1       ! MAXIMUM # OF ACTIVE SL+TW MODES IN ANY PHASE
#define NTWMMX 1       ! MAXIMUM # OF ACTIVE TWIN MODES IN ANY PHASE

#ifdef  TWO_SIGN_SLIP_SYSTEMS
  #define NSYSMX 24    ! MAXIMUM # OF ACTIVE SL+TW SYSTEMS IN ANY PHASE
#else
  #define NSYSMX 12    ! MAXIMUM # OF ACTIVE SL+TW SYSTEMS IN ANY PHASE
#endif


//===Max NR iterations allowed per voxel
//==Note: it may happen that for higher 
//==rate sensitivity params this goes to 
//=to even 1000 for some crystals if the 
//== LU decompoation is used and NOT 
//==quadraple precision. However, the 
//==GJE inverse implemented herein, takes
//==care of all! ;)
#define MAX_ITER_NR 512


// AbsoluteNoOutput results in the code not 
// printing any outputs to screen or file.
// Useful when profiling or scaling studies
#undef AbsoluteNoOutput
