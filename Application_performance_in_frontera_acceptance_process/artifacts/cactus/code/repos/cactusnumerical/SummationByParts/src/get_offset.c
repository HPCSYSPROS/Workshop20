#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

/***
 * A function to switch between CoordBase::GetBoundarySpecification and
 * MultiPatch::MultiPatch_GetBoundarySpecification, depending on whether
 * MultiPatch is present in compilation
 */
void get_shiftout ( const CCTK_POINTER_TO_CONST cctkGH_,
                    CCTK_INT *shiftout )
{
    cGH const * restrict const cctkGH = cctkGH_;

    DECLARE_CCTK_PARAMETERS

    int i;
    CCTK_INT nboundaryzones[6];
    CCTK_INT is_internal[6];
    CCTK_INT is_staggered[6];

    for (i=0;i<6;i++) shiftout[i] = 0;
    if (use_shiftout) {
        if (CCTK_IsFunctionAliased ("MultiPatch_GetBoundarySpecification")) {
            MultiPatch_GetBoundarySpecification (
                MultiPatch_GetMap(cctkGH_),
                6, 
                nboundaryzones, 
                is_internal, 
                is_staggered, 
                shiftout);
        } else if (CCTK_IsFunctionAliased ("GetBoundarySpecification")) {
            GetBoundarySpecification (6, 
                nboundaryzones, 
                is_internal, 
                is_staggered, 
                shiftout);
        } else {
            CCTK_WARN (CCTK_WARN_ABORT, 
            "Thorns providing GetBoundarySpecification function not found.");
        }
    }
}

/***
 * And it's fortran callable version...
 */
void CCTK_FCALL CCTK_FNAME(get_shiftout)( const CCTK_POINTER_TO_CONST cctkGH_,
                    CCTK_INT *shiftout )
{
  get_shiftout ( cctkGH_, shiftout );
}

/***
 * This function returns the effective values of local index ranges,
 * taking into account ghost zones and boundary_shiftout_* values.
 * If the use_shiftout=no, boundary offsets are set to zero.
 * Index ranges are set in C convention (starting from 0). I.e.,
 * to iterate over the grid in C, one might use
 *   for(i=imin[0];i<imax[0];i++) 
 *   for(j=imin[1];j<imax[1];j++)
 *   for(k=imin[2];k<imax[2];k++) 
 *       {...}
 * In Fortran, one would write the main loop like this:
 *   do i = imin(1)+1,imax(1)
 *     do j = imin(2)+1,imax(2)
 *       do k = imin(3)+1,imax(3)
 *       ...
 *       enddo
 *     enddo
 *   enddo
 * This function is declared in interface.ccl as GetLshIndexRanges 
 * to be used by other thorns, which need to know the actual index 
 * ranges for the region where SBP thorn will apply derivative /
 * dissipation operators.
 */
void get_lsh_iranges ( const CCTK_POINTER_TO_CONST cctkGH_, 
                       CCTK_INT *imin, 
                       CCTK_INT *imax)
{
    cGH const * restrict const cctkGH = cctkGH_;
    DECLARE_CCTK_PARAMETERS
    DECLARE_CCTK_ARGUMENTS

    int i;
    CCTK_INT shiftout[6];
    /*
    CCTK_INT nboundaryzones[6];
    CCTK_INT is_internal[6];
    CCTK_INT is_staggered[6];

    for (i=0;i<6;i++) shiftout[i] = 0;
    if (use_shiftout) {
        GetBoundarySpecification (6, 
            nboundaryzones, 
            is_internal, 
            is_staggered, 
            shiftout);
    }
    */
    get_shiftout (cctkGH_, shiftout);

    for (i=0;i<3;i++) {
        imin[i] = ((cctk_bbox[2*i]) ? shiftout[2*i] : cctk_nghostzones[i]);
        imax[i] = cctk_lsh[i] 
                - ((cctk_bbox[2*i+1]) ? shiftout[2*i+1] : cctk_nghostzones[i] );
    }

}


