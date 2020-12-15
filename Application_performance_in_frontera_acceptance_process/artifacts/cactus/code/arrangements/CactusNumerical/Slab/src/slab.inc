! -*-F90-*-

interface
   subroutine Slab_Transfer (ierr, cctkGH, dim, &
        src_gsh, src_lbnd, src_lsh, src_ash, &
        src_lbbox, src_ubbox, src_nghostzones, &
        src_off, src_str, src_len, &
        dst_gsh, dst_lbnd, dst_lsh, dst_ash, &
        dst_lbbox, dst_ubbox, dst_nghostzones, &
        dst_off, dst_str, dst_len, &
        xpose, flip, &
        options, &
        srctype, srcptr, &
        dsttype, dstptr)
     implicit none
     integer ierr
     CCTK_POINTER cctkGH 
     integer dim
     integer src_gsh(dim), src_lbnd(dim), src_lsh(dim), src_ash(dim)
     integer src_lbbox(dim), src_ubbox(dim), src_nghostzones(dim)
     integer src_off(dim), src_str(dim), src_len(dim)
     integer dst_gsh(dim), dst_lbnd(dim), dst_lsh(dim), dst_ash(dim)
     integer dst_lbbox(dim), dst_ubbox(dim), dst_nghostzones(dim)
     integer dst_off(dim), dst_str(dim), dst_len(dim)
     integer xpose(dim), flip(dim)
     integer options
     integer srctype
     CCTK_REAL srcptr(*)
     integer dsttype
     CCTK_REAL dstptr(*)
   end subroutine Slab_Transfer
end interface
