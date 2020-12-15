#define _inline_sse_mult_su3_mat_hwvec(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm0 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "movhps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((bb)->h[0].c[0]), \
                      "m" ((bb)->h[0].c[1]), \
                      "m" ((bb)->h[0].c[2]), \
                      "m" ((bb)->h[1].c[0]), \
                      "m" ((bb)->h[1].c[1]), \
                      "m" ((bb)->h[1].c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movss %0, %%xmm3 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm4 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "movss %4, %%xmm5 \n\t" \
                      "shufps $0x00, %%xmm3, %%xmm3 \n\t" \
                      "shufps $0x00, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm4, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "shufps $0x00, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm5, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %5, %%xmm6" \
                      : \
                      : \
                      "m" ((aa)->e[0][0].real), \
                      "m" ((aa)->e[0][1].real), \
                      "m" ((aa)->e[1][0].real), \
                      "m" ((aa)->e[1][2].real), \
                      "m" ((aa)->e[2][0].real), \
                      "m" ((aa)->e[2][1].real)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movss %0, %%xmm7 \n\t" \
                      "shufps $0x00, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm3 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x00, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "movss %3, %%xmm6 \n\t" \
                      "movss %4, %%xmm7 \n\t" \
                      "shufps $0xb1, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "shufps $0x00, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm7, %%xmm7 \n\t" \
                      "xorps %5, %%xmm0" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].real), \
                      "m" ((aa)->e[1][1].real), \
                      "m" ((aa)->e[2][2].real), \
                      "m" ((aa)->e[0][0].imag), \
                      "m" ((aa)->e[1][1].imag), \
                      "m" (_sse_sgn13)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("xorps %0, %%xmm1 \n\t" \
                      "xorps %1, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %2, %%xmm6 \n\t" \
                      "movss %3, %%xmm7 \n\t" \
                      "shufps $0x00, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      "movss %4, %%xmm6 \n\t" \
                      "movss %5, %%xmm7" \
                      : \
                      : \
                      "m" (_sse_sgn13), \
                      "m" (_sse_sgn13), \
                      "m" ((aa)->e[2][2].imag), \
                      "m" ((aa)->e[1][0].imag), \
                      "m" ((aa)->e[0][1].imag), \
                      "m" ((aa)->e[2][0].imag)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("shufps $0x00, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm3 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "movss %0, %%xmm0 \n\t" \
                      "movss %1, %%xmm6 \n\t" \
                      "movss %2, %%xmm7 \n\t" \
                      "shufps $0x00, %%xmm0, %%xmm0 \n\t" \
                      "shufps $0x00, %%xmm6, %%xmm6 \n\t" \
                      "shufps $0x00, %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm1, %%xmm6 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "addps %%xmm0, %%xmm3 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm4 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[0][2].imag), \
                      "m" ((aa)->e[2][1].imag), \
                      "m" ((aa)->e[1][2].imag)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      "movhps %%xmm3, %3 \n\t" \
                      "movhps %%xmm4, %4 \n\t" \
                      "movhps %%xmm5, %5" \
                      : \
                      "=m" ((cc)->h[0].c[0]), \
                      "=m" ((cc)->h[0].c[1]), \
                      "=m" ((cc)->h[0].c[2]), \
                      "=m" ((cc)->h[1].c[0]), \
                      "=m" ((cc)->h[1].c[1]), \
                      "=m" ((cc)->h[1].c[2])\
                      : \
                      : \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
}
