# hq_resonances
This is a repository where the FORTRAN code related to MCPR's work can be found.

Non-relativistic Hartree-Fock programs for heavy multiquark systems.
*   Discrete wavefunction as variational parameter
*   Equispaced grid. Simpson integration.
*   Charm and bottom masses fitted to Particle Data Group Experimental results

To run the code one must compile eigenrg.f previously. This is the code used for matrix diagonalization.

### File distribution:
*   GTO_basis directory: preliminary programs for He and $H_2$ systems. GTO basis functions are used. Hence, in these programs the variational parameter is not the wavefunction itself.
*   results directory: All results from .f90 programs are saved here.
*   HF_simp_* files: Main code for heavy quark systems.
    *   bc.f90: Mesons
    *   ccbb.f90: Differently flavoured tetraquarks
    *   cccc.f90: Same flavoured tetraquarks
    *   deuteron.f90: Preliminar program for the ccb-bbc hexaquark