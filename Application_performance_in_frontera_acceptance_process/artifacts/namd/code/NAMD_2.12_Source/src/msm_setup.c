/* msm_setup.c */

#include "msm_defn.h"


/* called by NL_msm_destroy() */
void NL_msm_cleanup(NL_Msm *pm) {
  int i;
#ifdef NL_USE_CUDA
  if (pm->msmflags & NL_MSM_COMPUTE_CUDA_GRID_CUTOFF) {
    NL_msm_cuda_cleanup_gridcutoff(pm);
  }
#endif /* NL_USE_CUDA */
  if (pm->msmflags & NL_MSM_COMPUTE_SPREC) {
    for (i = 0;  i < pm->maxlevels;  i++) {
      GRID_DONE( &(pm->qh_f[i]) );
      GRID_DONE( &(pm->eh_f[i]) );
      GRID_DONE( &(pm->gc_f[i]) );
    }
  }
  else {
    for (i = 0;  i < pm->maxlevels;  i++) {
      GRID_DONE( &(pm->qh[i]) );
      GRID_DONE( &(pm->eh[i]) );
      GRID_DONE( &(pm->gc[i]) );
    }
  }
  free(pm->lzd);
  free(pm->lyzd);
  free(pm->lzd_f);
  free(pm->lyzd_f);
  free(pm->qh);
  free(pm->eh);
  free(pm->gc);
  free(pm->qh_f);
  free(pm->eh_f);
  free(pm->gc_f);
}


static int setup_bins_1away(NL_Msm *pm);
static int setup_bins_k_away(NL_Msm *pm);

static int setup_cell_vectors(NL_Msm *pm);

static int setup_grids(NL_Msm *pm);

/* called by setup_grids() */
static int setup_hgrid_1d(
    NL_Msm *pm,
    double len,           /* cell length */
    double *hh,           /* determine h */
    int *nn,              /* determine number grid spacings covering cell */
    int *aindex,          /* determine smallest grid index */
    int *bindex,          /* determine largest grid index */
    int isperiodic        /* is this dimension periodic? */
    );

static int print_status(NL_Msm *msm);


int NL_msm_setup(
    NL_Msm *pm,
    double cutoff,
    double cellvec1[3],
    double cellvec2[3],
    double cellvec3[3],
    double cellcenter[3],
    int msmflags
    ) {

  int rc = 0;

  /* check sanity of msmflags */
  if ((~NL_MSM_ALL_FLAGS & msmflags) != 0 ||
      (NL_MSM_COMPUTE_ALL & msmflags) == 0) {
    return NL_MSM_ERROR_PARAM;
  }

  /* report timings? */
  pm->report_timings = ((msmflags & NL_MSM_REPORT_TIMINGS) != 0);

  /* for now, permit only orthogonal cell aligned with coordinate axes */
  if (cellvec1[1] != 0 || cellvec1[2] != 0 ||
      cellvec2[0] != 0 || cellvec2[2] != 0 ||
      cellvec3[0] != 0 || cellvec3[1] != 0 ||
      cellvec1[0] <= 0 || cellvec2[1] <= 0 || cellvec3[2] <= 0) {
    return NL_MSM_ERROR_SUPPORT;
  }

  /* check sanity of cutoff wrt expected cell;
   * XXX cell widths must be at least the cutoff distance length */
  if (cutoff <= 0 ||
      (cutoff > cellvec1[0] && (msmflags & NL_MSM_PERIODIC_VEC1)) ||
      (cutoff > cellvec2[1] && (msmflags & NL_MSM_PERIODIC_VEC2)) ||
      (cutoff > cellvec3[2] && (msmflags & NL_MSM_PERIODIC_VEC3))) {
    return NL_MSM_ERROR_PARAM;
  }

  pm->msmflags = msmflags;
  pm->cellvec1[0] = cellvec1[0];
  pm->cellvec1[1] = cellvec1[1];
  pm->cellvec1[2] = cellvec1[2];
  pm->cellvec2[0] = cellvec2[0];
  pm->cellvec2[1] = cellvec2[1];
  pm->cellvec2[2] = cellvec2[2];
  pm->cellvec3[0] = cellvec3[0];
  pm->cellvec3[1] = cellvec3[1];
  pm->cellvec3[2] = cellvec3[2];
  pm->cellcenter[0] = cellcenter[0];
  pm->cellcenter[1] = cellcenter[1];
  pm->cellcenter[2] = cellcenter[2];

  pm->a = cutoff;

  rc = setup_cell_vectors(pm);
  if (rc) return rc;

  if (msmflags & NL_MSM_COMPUTE_SHORT_RANGE) {
    /* set up bins for short-range part */
    if (msmflags & NL_MSM_COMPUTE_1AWAY) {
      rc = setup_bins_1away(pm);
    }
    else {
      rc = setup_bins_k_away(pm);
    }
    if (rc) return rc;
  }

  if (msmflags & NL_MSM_COMPUTE_LONG_RANGE) {
    /* set up grid hierarchy for long-range part */
    rc = setup_grids(pm);
    if (rc) return rc;
  }

  if (msmflags & NL_MSM_COMPUTE_SPREC) {
    /* fill out single precision data if needed */
    pm->cellvec1_f[0] = (float) pm->cellvec1[0];
    pm->cellvec1_f[1] = (float) pm->cellvec1[1];
    pm->cellvec1_f[2] = (float) pm->cellvec1[2];
    pm->cellvec2_f[0] = (float) pm->cellvec2[0];
    pm->cellvec2_f[1] = (float) pm->cellvec2[1];
    pm->cellvec2_f[2] = (float) pm->cellvec2[2];
    pm->cellvec3_f[0] = (float) pm->cellvec3[0];
    pm->cellvec3_f[1] = (float) pm->cellvec3[1];
    pm->cellvec3_f[2] = (float) pm->cellvec3[2];
    pm->cellcenter_f[0] = (float) pm->cellcenter[0];
    pm->cellcenter_f[1] = (float) pm->cellcenter[1];
    pm->cellcenter_f[2] = (float) pm->cellcenter[2];
    pm->recipvec1_f[0] = (float) pm->recipvec1[0];
    pm->recipvec1_f[1] = (float) pm->recipvec1[1];
    pm->recipvec1_f[2] = (float) pm->recipvec1[2];
    pm->recipvec2_f[0] = (float) pm->recipvec2[0];
    pm->recipvec2_f[1] = (float) pm->recipvec2[1];
    pm->recipvec2_f[2] = (float) pm->recipvec2[2];
    pm->recipvec3_f[0] = (float) pm->recipvec3[0];
    pm->recipvec3_f[1] = (float) pm->recipvec3[1];
    pm->recipvec3_f[2] = (float) pm->recipvec3[2];
    pm->hx_f = pm->hx;
    pm->hy_f = pm->hy;
    pm->hz_f = pm->hz;
    pm->a_f = pm->a;
    pm->gx_f = pm->gx;
    pm->gy_f = pm->gy;
    pm->gz_f = pm->gz;
  }

  print_status(pm);

#ifdef NL_USE_CUDA
  if (msmflags & NL_MSM_COMPUTE_CUDA_GRID_CUTOFF) {
    rc = NL_msm_cuda_setup_gridcutoff(pm);
    if (rc == NL_MSM_SUCCESS) {
      printf("Using CUDA for grid cutoff computation\n");
    }
    else {
      printf("Unable to set up CUDA for grid cutoff computation\n");
      if (msmflags & NL_MSM_COMPUTE_CUDA_FALL_BACK) {
        NL_msm_cuda_cleanup_gridcutoff(pm);
        printf("Falling back on CPU\n");
       	pm->msmflags &= ~NL_MSM_COMPUTE_CUDA_GRID_CUTOFF;
      }
      else return rc;
    }
  }
#else
  if (msmflags & NL_MSM_COMPUTE_CUDA_GRID_CUTOFF) {
    if (msmflags & NL_MSM_COMPUTE_CUDA_FALL_BACK) {
      printf("Falling back on CPU\n");
      pm->msmflags &= ~NL_MSM_COMPUTE_CUDA_GRID_CUTOFF;
    }
    else return NL_MSM_ERROR_SUPPORT;
  }
#endif /* NL_USE_CUDA */

  return NL_MSM_SUCCESS;
}


typedef struct InterpParams_t {
  int nu;
  int stencil;
  int omega;
} InterpParams;

static InterpParams INTERP_PARAMS[] = {
  { 1, 4, 6 },    /* cubic */
  { 2, 6, 10 },   /* quintic */
  { 2, 6, 10 },   /* quintic, C2 */
  { 3, 8, 14 },   /* septic */
  { 3, 8, 14 },   /* septic, C3 */
  { 4, 10, 18 },  /* nonic */
  { 4, 10, 18 },  /* nonic, C4 */
  { 1, 4, 6 },    /* B-spline */
};


int print_status(NL_Msm *pm) {
  int k;
  int ispx = (pm->msmflags & NL_MSM_PERIODIC_VEC1);
  int ispy = (pm->msmflags & NL_MSM_PERIODIC_VEC2);
  int ispz = (pm->msmflags & NL_MSM_PERIODIC_VEC3);
  int ispany = (pm->msmflags & NL_MSM_PERIODIC_ALL);
  int ispall = (ispany == NL_MSM_PERIODIC_ALL);

  const double xlen = pm->cellvec1[0];  /* XXX */
  const double ylen = pm->cellvec2[1];
  const double zlen = pm->cellvec3[2];

  printf("#MSM SETUP\n");
  printf("#  cell lengths=  %g  %g  %g\n", xlen, ylen, zlen);
  printf("#  grid origin=  %g  %g  %g\n", pm->gx, pm->gy, pm->gz);
  if (ispall) {
    printf("#  periodic boundaries\n");
  }
  else if (!ispany) {
    printf("#  non-periodic boundaries\n");
  }
  else {
    printf("#  periodic boundaries in:%s%s%s\n",
        (ispx ? "  x" : ""), (ispy ? "  y" : ""), (ispz ? "  z" : ""));
  }
  printf("#  cutoff= %g\n", pm->a);
  printf("#  grid spacing= %g\n", pm->gridspacing);
  printf("#  hx, hy, hz= %g, %g, %g\n", pm->hx, pm->hy, pm->hz);
  printf("#  h-grid size= %d x %d x %d\n", pm->nx, pm->ny, pm->nz);
  printf("#  number of levels= %d\n", pm->nlevels);
  printf("#  approximation= %s\n", NL_msm_approx_name(pm->approx));
  printf("#  splitting= %s\n", NL_msm_split_name(pm->split));
  printf("#  grid hierarchy:\n");
  if (pm->msmflags & NL_MSM_COMPUTE_SPREC) {
    for (k = 0;  k < pm->nlevels;  k++) {
      NL_Msmgrid_float *qh = &(pm->qh_f[k]);
      int ia = qh->i0;
      int ib = ia + qh->ni - 1;
      int ja = qh->j0;
      int jb = ja + qh->nj - 1;
      int ka = qh->k0;
      int kb = ka + qh->nk - 1;
      printf("#  level= %d:  [%d..%d] x [%d..%d] x [%d..%d]\n",
          k, ia, ib, ja, jb, ka, kb);
    }
  }
  else {
    for (k = 0;  k < pm->nlevels;  k++) {
      NL_Msmgrid_double *qh = &(pm->qh[k]);
      int ia = qh->i0;
      int ib = ia + qh->ni - 1;
      int ja = qh->j0;
      int jb = ja + qh->nj - 1;
      int ka = qh->k0;
      int kb = ka + qh->nk - 1;
      printf("#  level= %d:  [%d..%d] x [%d..%d] x [%d..%d]\n",
          k, ia, ib, ja, jb, ka, kb);
    }
  }
  return NL_MSM_SUCCESS;
}


int setup_cell_vectors(NL_Msm *pm) {
  double *u = pm->cellvec1;
  double *v = pm->cellvec2;
  double *w = pm->cellvec3;
  double *bu = pm->recipvec1;
  double *bv = pm->recipvec2;
  double *bw = pm->recipvec3;
  double c[3], s;

  c[X] = v[Y]*w[Z] - v[Z]*w[Y];  /* v CROSS w */
  c[Y] = v[Z]*w[X] - v[X]*w[Z];
  c[Z] = v[X]*w[Y] - v[Y]*w[X];
  s = 1 / (u[X]*c[X] + u[Y]*c[Y] + u[Z]*c[Z]);  /* 1 / (c DOT u) */
  bu[X] = s*c[X];
  bu[Y] = s*c[Y];
  bu[Z] = s*c[Z];

  c[X] = w[Y]*u[Z] - w[Z]*u[Y];  /* w CROSS u */
  c[Y] = w[Z]*u[X] - w[X]*u[Z];
  c[Z] = w[X]*u[Y] - w[Y]*u[X];
  s = 1 / (v[X]*c[X] + v[Y]*c[Y] + v[Z]*c[Z]);  /* 1 / (c DOT v) */
  bv[X] = s*c[X];
  bv[Y] = s*c[Y];
  bv[Z] = s*c[Z];

  c[X] = u[Y]*v[Z] - u[Z]*v[Y];  /* u CROSS v */
  c[Y] = u[Z]*v[X] - u[X]*v[Z];
  c[Z] = u[X]*v[Y] - u[Y]*v[X];
  s = 1 / (w[X]*c[X] + w[Y]*c[Y] + w[Z]*c[Z]);  /* 1 / (c DOT w) */
  bw[X] = s*c[X];
  bw[Y] = s*c[Y];
  bw[Z] = s*c[Z];

  return NL_MSM_SUCCESS;
}


int setup_bins_1away(NL_Msm *pm) {
  const double *u = pm->cellvec1;
  const double *v = pm->cellvec2;
  const double *w = pm->cellvec3;
  double *ru = pm->recipvec1;
  double *rv = pm->recipvec2;
  double *rw = pm->recipvec3;
  double p[3];
  double pulen, pvlen, pwlen, s;
  double cutoff = pm->a;
  int nbx, nby, nbz, numbins;
  int ispx = ((pm->msmflags & NL_MSM_PERIODIC_VEC1) != 0);
  int ispy = ((pm->msmflags & NL_MSM_PERIODIC_VEC2) != 0);
  int ispz = ((pm->msmflags & NL_MSM_PERIODIC_VEC3) != 0);

  /* calculate the number of atom bins in each basis vector dimension,
   * such that each bin (a parallelepiped) inscribes the cutoff-cube;
   * for periodic boundaries, we have to choose equal sized bins with
   * length >= cutoff that tile the cell; for non-periodic boundaries,
   * we can have bins of length = cutoff */

  /* find the largest orthogonal box inscribed within parallelepiped cell
   * by taking orthogonal projections onto cross products of basis vectors */

  p[X] = v[Y]*w[Z] - v[Z]*w[Y];  /* p = v CROSS w */
  p[Y] = v[Z]*w[X] - v[X]*w[Z];
  p[Z] = v[X]*w[Y] - v[Y]*w[X];
  /* s = (u DOT p) / (p DOT p) */
  s = (u[X]*p[X] + u[Y]*p[Y] + u[Z]*p[Z]) / (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  p[X] *= s;  /* p is orthogonal projection of u onto v CROSS w */
  p[Y] *= s;
  p[Z] *= s;
  pulen = sqrt(p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);

  p[X] = w[Y]*u[Z] - w[Z]*u[Y];  /* p = w CROSS u */
  p[Y] = w[Z]*u[X] - w[X]*u[Z];
  p[Z] = w[X]*u[Y] - w[Y]*u[X];
  /* s = (v DOT p) / (p DOT p) */
  s = (v[X]*p[X] + v[Y]*p[Y] + v[Z]*p[Z]) / (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  p[X] *= s;  /* p is orthogonal projection of v onto w CROSS u */
  p[Y] *= s;
  p[Z] *= s;
  pvlen = sqrt(p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);

  p[X] = u[Y]*v[Z] - u[Z]*v[Y];  /* p = u CROSS v */
  p[Y] = u[Z]*v[X] - u[X]*v[Z];
  p[Z] = u[X]*v[Y] - u[Y]*v[X];
  /* s = (w DOT p) / (p DOT p) */
  s = (w[X]*p[X] + w[Y]*p[Y] + w[Z]*p[Z]) / (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  p[X] *= s;  /* p is orthogonal projection of w onto u CROSS v */
  p[Y] *= s;
  p[Z] *= s;
  pwlen = sqrt(p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);

  /* find the number of cells along each dimension by dividing the
   * orthogonal projection vector lengths by the cutoff distance */
  nbx = (ispx ? (int) floor(pulen / cutoff) : (int) ceil(pulen / cutoff));
  nby = (ispy ? (int) floor(pvlen / cutoff) : (int) ceil(pvlen / cutoff));
  nbz = (ispz ? (int) floor(pwlen / cutoff) : (int) ceil(pwlen / cutoff));
  if (nbx == 0 || nby == 0 || nbz == 0) {
    return NL_MSM_ERROR_PARAM;  /* cutoff is too large for cell basis */
  }
  numbins = nbx * nby * nbz;
  /* we don't know the number of atoms until compute */

  /* allocate one atom index for each bin */
  if (pm->maxbins < numbins) {
    void *vptr = malloc(numbins * sizeof(int));
    if (vptr == NULL) return NL_MSM_ERROR_MALLOC;
    pm->bin = (int *) vptr;
    pm->maxbins = numbins;
  }
  pm->nbx = nbx;
  pm->nby = nby;
  pm->nbz = nbz;
  pm->numbins = numbins;

  /* must allocate next index array when we know number of atoms */

  /* rescale recip space vectors for non-periodic expansion */
  if (ispx == 0) {
    s = pulen / (nbx * cutoff);
    ru[X] *= s;
    ru[Y] *= s;
    ru[Z] *= s;
  }
  if (ispy == 0) {
    s = pvlen / (nby * cutoff);
    rv[X] *= s;
    rv[Y] *= s;
    rv[Z] *= s;
  }
  if (ispz == 0) {
    s = pwlen / (nbz * cutoff);
    rw[X] *= s;
    rw[Y] *= s;
    rw[Z] *= s;
  }
  return NL_MSM_SUCCESS;
}


int setup_bins_k_away(NL_Msm *pm) {
  double *u = pm->cellvec1;
  double *v = pm->cellvec2;
  double *w = pm->cellvec3;
  double *ru = pm->recipvec1;
  double *rv = pm->recipvec2;
  double *rw = pm->recipvec3;
  double *bu = pm->bu;
  double *bv = pm->bv;
  double *bw = pm->bw;
  double p[3];
  double pulen, pvlen, pwlen, s;
  double cutoff = pm->a;
  int nbx, nby, nbz, numbins;
  int ispx = ((pm->msmflags & NL_MSM_PERIODIC_VEC1) != 0);
  int ispy = ((pm->msmflags & NL_MSM_PERIODIC_VEC2) != 0);
  int ispz = ((pm->msmflags & NL_MSM_PERIODIC_VEC3) != 0);
  int nbrhdmax, nradius, ndiameter, index, i, j, k;
  double volume;
  double binvol, abincnt, binlen, c;
  double min2;

  /* find the largest orthogonal box inscribed within parallelepiped cell
   * by taking orthogonal projections onto cross products of basis vectors */
  p[X] = v[Y]*w[Z] - v[Z]*w[Y];  /* p = v CROSS w */
  p[Y] = v[Z]*w[X] - v[X]*w[Z];
  p[Z] = v[X]*w[Y] - v[Y]*w[X];
  /* s = (u DOT p) / (p DOT p) */
  s = (u[X]*p[X] + u[Y]*p[Y] + u[Z]*p[Z]) / (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  p[X] *= s;  /* p is orthogonal projection of u onto v CROSS w */
  p[Y] *= s;
  p[Z] *= s;
  pulen = sqrt(p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);

  p[X] = w[Y]*u[Z] - w[Z]*u[Y];  /* p = w CROSS u */
  p[Y] = w[Z]*u[X] - w[X]*u[Z];
  p[Z] = w[X]*u[Y] - w[Y]*u[X];
  /* s = (v DOT p) / (p DOT p) */
  s = (v[X]*p[X] + v[Y]*p[Y] + v[Z]*p[Z]) / (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  p[X] *= s;  /* p is orthogonal projection of v onto w CROSS u */
  p[Y] *= s;
  p[Z] *= s;
  pvlen = sqrt(p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);

  p[X] = u[Y]*v[Z] - u[Z]*v[Y];  /* p = u CROSS v */
  p[Y] = u[Z]*v[X] - u[X]*v[Z];
  p[Z] = u[X]*v[Y] - u[Y]*v[X];
  volume = w[X]*p[X] + w[Y]*p[Y] + w[Z]*p[Z];
  /* s = (w DOT p) / (p DOT p) */
  s = volume / (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  p[X] *= s;  /* p is orthogonal projection of w onto u CROSS v */
  p[Y] *= s;
  p[Z] *= s;
  pwlen = sqrt(p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);

  if ((ispx && cutoff > pulen) ||
      (ispy && cutoff > pvlen) ||
      (ispz && cutoff > pwlen)) {
    printf("Cutoff %.3g too big for cell along basis%s%s%s\n",
        cutoff, (cutoff > pulen ? " x" : ""), (cutoff > pvlen ?  " y" : ""),
        (cutoff > pwlen ?  " z" : ""));
    return NL_MSM_ERROR_PARAM;
  }

  /* calculate the ideal bin volume based on a fixed bin size (nbinslots),
   * the particle density (density), and a desired fill ratio (binfill) */
  binvol = pm->binfill * pm->nbinslots / pm->density;
  abincnt = volume / binvol;
  binlen = pow(binvol, 1./3);  /* ideal side length - use for nonperiodic */

  /* factor "abincnt" into 3 parts, each part proportional to the
   * lengths of the orthogonal projection vectors calculated above */
  c = pow(abincnt / (pulen*pvlen*pwlen), 1./3);  /* proportionality const */
  nbx = (int) ceil(c * pulen);
  nby = (int) ceil(c * pvlen);
  nbz = (int) ceil(c * pwlen);
  numbins = nbx * nby * nbz;

  printf("nbx=%d  nby=%d  nbz=%d  numbins=%d\n", nbx, nby, nbz, numbins);

  /* allocate one atom index for each bin */
  if (pm->maxbins < numbins) {
    void *vptr = malloc(numbins * sizeof(int));
    if (vptr == NULL) return NL_MSM_ERROR_MALLOC;
    pm->bin = (int *) vptr;
    pm->maxbins = numbins;
  }
  pm->nbx = nbx;
  pm->nby = nby;
  pm->nbz = nbz;
  pm->numbins = numbins;

  /* rescale basis and recip space vectors for non-periodic expansion */
  if ( ! ispx) {
    s = pulen / (nbx * binlen);
    ru[X] *= s;
    ru[Y] *= s;
    ru[Z] *= s;
    s = 1./s;
    u[X] *= s;
    u[Y] *= s;
    u[Z] *= s;
  }
  if ( ! ispy) {
    s = pvlen / (nby * binlen);
    rv[X] *= s;
    rv[Y] *= s;
    rv[Z] *= s;
    s = 1./s;
    v[X] *= s;
    v[Y] *= s;
    v[Z] *= s;
  }
  if ( ! ispz) {
    s = pwlen / (nbz * binlen);
    rw[X] *= s;
    rw[Y] *= s;
    rw[Z] *= s;
    s = 1./s;
    w[X] *= s;
    w[Y] *= s;
    w[Z] *= s;
  }

  /* scale basis vectors by number of bins to get bin basis vectors */
  s = 1./nbx;
  bu[X] = s * u[X];
  bu[Y] = s * u[Y];
  bu[Z] = s * u[Z];
  s = 1./nby;
  bv[X] = s * v[X];
  bv[Y] = s * v[Y];
  bv[Z] = s * v[Z];
  s = 1./nbz;
  bw[X] = s * w[X];
  bw[Y] = s * w[Y];
  bw[Z] = s * w[Z];

  /* determine the neighborhood of bins */

  /* first find minimum width of cell */
  min2 = (bu[X]*bu[X] + bu[Y]*bu[Y] + bu[Z]*bu[Z]);
  s = (bv[X]*bv[X] + bv[Y]*bv[Y] + bv[Z]*bv[Z]);
  if (min2 > s) min2 = s;
  s = (bw[X]*bw[X] + bw[Y]*bw[Y] + bw[Z]*bw[Z]);
  if (min2 > s) min2 = s;

  /* also find the minimum of the four major diagonals */
  p[X] = bu[X] + bv[X] + bw[X];
  p[Y] = bu[Y] + bv[Y] + bw[Y];
  p[Z] = bu[Z] + bv[Z] + bw[Z];
  s = (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  if (min2 > s) min2 = s;
  p[X] = -bu[X] + bv[X] + bw[X];
  p[Y] = -bu[Y] + bv[Y] + bw[Y];
  p[Z] = -bu[Z] + bv[Z] + bw[Z];
  s = (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  if (min2 > s) min2 = s;
  p[X] = bu[X] - bv[X] + bw[X];
  p[Y] = bu[Y] - bv[Y] + bw[Y];
  p[Z] = bu[Z] - bv[Z] + bw[Z];
  s = (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  if (min2 > s) min2 = s;
  p[X] = bu[X] + bv[X] - bw[X];
  p[Y] = bu[Y] + bv[Y] - bw[Y];
  p[Z] = bu[Z] + bv[Z] - bw[Z];
  s = (p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
  if (min2 > s) min2 = s;

  /* the neighborhood is contained in the cube of nradius bins,
   * store only upper half of bin neighborhood */
  nradius = (int) ceil(sqrt(cutoff*cutoff / min2));
  ndiameter = 2 * nradius + 1;
  nbrhdmax = (ndiameter * ndiameter * ndiameter) / 2 + 1;

#if 0
  printf("Neighborhood radius = %d\n", nradius);
  printf("min2 = %g\n", min2);
#endif

  /* allocate space for the entire cube */
  if (pm->nbrhd) free(pm->nbrhd);
  pm->nbrhd = (int *) calloc(3*nbrhdmax, sizeof(int));
  if (pm->nbrhd == NULL) return NL_MSM_ERROR_MALLOC;
  pm->nbrhdmax = nbrhdmax;
  pm->nradius = nradius;  /* save neighborhood radius for diagnostic purposes */

  /* trim the neighborhood */
  index = 0;
  for (k = -nradius;  k <= nradius;  k++) {
    for (j = -nradius;  j <= nradius;  j++) {
      for (i = -nradius;  i <= nradius;  i++) {
        /* XXX should we remove (0,0,0) from bin neighborhood? */
        int ip, jp, kp, iq, jq, kq;
        int is_within = 0;
        int binindex = (k * ndiameter + j) * ndiameter + i;
        if (binindex < 0) continue;
        for (kp = 0;  kp <= 1;  kp++) {
          for (jp = 0;  jp <= 1;  jp++) {
            for (ip = 0;  ip <= 1;  ip++) {
              p[X] = (i+ip)*bu[X] + (j+jp)*bv[X] + (k+kp)*bw[X];
              p[Y] = (i+ip)*bu[Y] + (j+jp)*bv[Y] + (k+kp)*bw[Y];
              p[Z] = (i+ip)*bu[Z] + (j+jp)*bv[Z] + (k+kp)*bw[Z];
              for (kq = 0;  kq <= 1;  kq++) {
                for (jq = 0;  jq <= 1;  jq++) {
                  for (iq = 0;  iq <= 1;  iq++) {
                    double q[3];
                    double dx, dy, dz, r2;
                    q[X] = iq*bu[X] + jq*bv[X] + kq*bw[X];
                    q[Y] = iq*bu[Y] + jq*bv[Y] + kq*bw[Y];
                    q[Z] = iq*bu[Z] + jq*bv[Z] + kq*bw[Z];
                    dx = p[X] - q[X];
                    dy = p[Y] - q[Y];
                    dz = p[Z] - q[Z];
                    r2 = dx*dx + dy*dy + dz*dz;
                    is_within |= (r2 < cutoff*cutoff);
                  }
                }
              } /* end q loop */

            }
          }
        } /* end p loop */
        if (is_within) {
          pm->nbrhd[index++] = i;
          pm->nbrhd[index++] = j;
          pm->nbrhd[index++] = k;
        }

      }
    }
  }
  pm->nbrhdlen = index / 3;

  return NL_MSM_SUCCESS;
}


int setup_hgrid_1d(
    NL_Msm *pm,
    double len,           /* cell length */
    double *hh,           /* determine h */
    int *nn,              /* determine number grid spacings covering cell */
    int *aindex,          /* determine smallest grid index */
    int *bindex,          /* determine largest grid index */
    int isperiodic        /* is this dimension periodic? */
    ) {

  const int nu = INTERP_PARAMS[pm->approx].nu;  /* interp stencil radius */

  ASSERT(hmax > 0);
  if (isperiodic) {
    const double hmin = (4./5) * pm->gridspacing;  /* minimum bound on h */
    const double hmax = 1.5 * hmin;
    double h = len;
    int n = 1;    /* start with one grid point across domain */
    while (h >= hmax) {
      h *= 0.5;   /* halve h */
      n <<= 1;    /* double grid points */
    }
    if (h < hmin) {
      if (n < 4) {  /* either len is too small or hmin is too large */
        return NL_MSM_ERROR_PARAM;
      }
      h *= (4./3);  /* scale h by 4/3 */
      n >>= 2;      /* scale n by 3/4 */
      n *= 3;
    }
    /* now we have:  hmin <= h < hmax */
    /* now we have:  n is power of two times no more than one power of 3 */
    *hh = h;
    *nn = n;
    *aindex = 0;
    *bindex = n-1;
  }
  else {  /* non-periodic */
    double h = pm->gridspacing;
    int n = (int) floorf(len / h) + 1;
    *hh = h;
    *nn = n;
    *aindex = -nu;
    *bindex = n + nu;
  }
  return NL_MSM_SUCCESS;
}


int setup_grids(NL_Msm *pm) {
  const int nu = INTERP_PARAMS[pm->approx].nu;
  const int omega = INTERP_PARAMS[pm->approx].omega;
  const int split = pm->split;
  const int ispx = (pm->msmflags & NL_MSM_PERIODIC_VEC1);
  const int ispy = (pm->msmflags & NL_MSM_PERIODIC_VEC2);
  const int ispz = (pm->msmflags & NL_MSM_PERIODIC_VEC3);
  const int ispany = (pm->msmflags & NL_MSM_PERIODIC_ALL);

  const int issprec = (pm->msmflags & NL_MSM_COMPUTE_SPREC);

  const double xlen = pm->cellvec1[0];  /* XXX */
  const double ylen = pm->cellvec2[1];
  const double zlen = pm->cellvec3[2];

  const double a = pm->a;
  double hx, hy, hz;
  double scaling;
  double d;  /* temporary for SPOLY derivative */

  NL_Msmgrid_double *p = NULL;
  NL_Msmgrid_float *p_f = NULL;
  int ia, ib, ja, jb, ka, kb, ni, nj, nk;
  int nx, ny, nz;  /* counts the grid points that span just the domain */

  int i, j, k, n;
  int index;
  int level, toplevel, nlevels, maxlevels;
  int lastnelems = 1;
  int isclamped = 0;
  int done, alldone;

  int rc = 0;

  rc = setup_hgrid_1d(pm, xlen, &hx, &nx, &ia, &ib, ispx);
  if (rc) return rc;

  rc = setup_hgrid_1d(pm, ylen, &hy, &ny, &ja, &jb, ispy);
  if (rc) return rc;

  rc = setup_hgrid_1d(pm, zlen, &hz, &nz, &ka, &kb, ispz);
  if (rc) return rc;

  pm->hx = hx;
  pm->hy = hy;
  pm->hz = hz;

  /* XXX set coordinate for h-grid (0,0,0) point */
  pm->gx = pm->cellcenter[0] - ((nx >> 1) * hx);
  pm->gy = pm->cellcenter[1] - ((ny >> 1) * hy);
  pm->gz = pm->cellcenter[2] - ((nz >> 1) * hz);

  pm->nx = nx;
  pm->ny = ny;
  pm->nz = nz;

  ni = ib - ia + 1;
  nj = jb - ja + 1;
  nk = kb - ka + 1;

  /* allocate temp buffer space for factored grid transfer */
  n = (nk > omega ? nk : omega);  /* row along z-dimension */
  if (pm->max_lzd < n) {
    if (issprec) {
      float *t;
      t = (float *) realloc(pm->lzd, n * sizeof(float));
      if (NULL == t) return NL_MSM_ERROR_MALLOC;
      pm->lzd_f = t;
    }
    else {
      double *t;
      t = (double *) realloc(pm->lzd, n * sizeof(double));
      if (NULL == t) return NL_MSM_ERROR_MALLOC;
      pm->lzd = t;
    }
    pm->max_lzd = n;
  }
  n *= (nj > omega ? nj : omega);  /* plane along yz-dimensions */
  if (pm->max_lyzd < n) {
    if (issprec) {
      float *t;
      t = (float *) realloc(pm->lyzd, n * sizeof(float));
      if (NULL == t) return NL_MSM_ERROR_MALLOC;
      pm->lyzd_f = t;
    }
    else {
      double *t;
      t = (double *) realloc(pm->lyzd, n * sizeof(double));
      if (NULL == t) return NL_MSM_ERROR_MALLOC;
      pm->lyzd = t;
    }
    pm->max_lyzd = n;
  }

  nlevels = pm->nlevels;
  if (nlevels <= 0) {
    /* automatically set number of levels */
    n = ni;
    if (n < nj) n = nj;
    if (n < nk) n = nk;
    for (maxlevels = 1;  n > 0;  n >>= 1)  maxlevels++;
    nlevels = maxlevels;
    if (ispany == 0) {  /* no periodicity */
      int omega3 = omega * omega * omega;
      int nhalf = (int) sqrtf(ni*nj*nk);  /* scale down for performance? */
      lastnelems = (nhalf > omega3 ? nhalf : omega3);
      isclamped = 1;
    }
  }
  else {
    /* user-defined number of levels */
    maxlevels = nlevels;
  }

  /* allocate any additional levels that may be needed */
  if (pm->maxlevels < maxlevels) {
    void *vqh, *veh, *vgc;
    if (issprec) {
      vqh = realloc(pm->qh, maxlevels * sizeof(NL_Msmgrid_float));
      if (NULL == vqh) return NL_MSM_ERROR_MALLOC;
      veh = realloc(pm->eh, maxlevels * sizeof(NL_Msmgrid_float));
      if (NULL == veh) return NL_MSM_ERROR_MALLOC;
      vgc = realloc(pm->gc, maxlevels * sizeof(NL_Msmgrid_float));
      if (NULL == vgc) return NL_MSM_ERROR_MALLOC;
      pm->qh_f = (NL_Msmgrid_float *) vqh;
      pm->eh_f = (NL_Msmgrid_float *) veh;
      pm->gc_f = (NL_Msmgrid_float *) vgc;
      /* initialize the newest grids appended to array */
      for (level = pm->maxlevels;  level < maxlevels;  level++) {
        GRID_INIT( &(pm->qh_f[level]) );
        GRID_INIT( &(pm->eh_f[level]) );
        GRID_INIT( &(pm->gc_f[level]) );
      }
    }
    else {
      vqh = realloc(pm->qh, maxlevels * sizeof(NL_Msmgrid_double));
      if (NULL == vqh) return NL_MSM_ERROR_MALLOC;
      veh = realloc(pm->eh, maxlevels * sizeof(NL_Msmgrid_double));
      if (NULL == veh) return NL_MSM_ERROR_MALLOC;
      vgc = realloc(pm->gc, maxlevels * sizeof(NL_Msmgrid_double));
      if (NULL == vgc) return NL_MSM_ERROR_MALLOC;
      pm->qh = (NL_Msmgrid_double *) vqh;
      pm->eh = (NL_Msmgrid_double *) veh;
      pm->gc = (NL_Msmgrid_double *) vgc;
      /* initialize the newest grids appended to array */
      for (level = pm->maxlevels;  level < maxlevels;  level++) {
        GRID_INIT( &(pm->qh[level]) );
        GRID_INIT( &(pm->eh[level]) );
        GRID_INIT( &(pm->gc[level]) );
      }
    }
    pm->maxlevels = maxlevels;
  }

  level = 0;
  done = 0;
  alldone = 0;
  do {
    if (issprec) {
      GRID_RESIZE( &(pm->qh_f[level]), float, ia, ni, ja, nj, ka, nk);
      GRID_RESIZE( &(pm->eh_f[level]), float, ia, ni, ja, nj, ka, nk);
    }
    else {
      GRID_RESIZE( &(pm->qh[level]), double, ia, ni, ja, nj, ka, nk);
      GRID_RESIZE( &(pm->eh[level]), double, ia, ni, ja, nj, ka, nk);
    }

    if (++level == nlevels)    done |= 0x07;  /* user limit on levels */

    alldone = (done == 0x07);  /* make sure all dimensions are done */

    if (isclamped) {
      int nelems = ni * nj * nk;
      if (nelems <= lastnelems)  done |= 0x07;
    }

    if (ispx) {
      ni >>= 1;
      ib = ni-1;
      if (ni & 1)              done |= 0x07;  /* == 3 or 1 */
      else if (ni == 2)        done |= 0x01;  /* can do one more */
    }
    else {
      ia = -((-ia+1)/2) - nu;
      ib = (ib+1)/2 + nu;
      ni = ib - ia + 1;
      if (ni <= omega)         done |= 0x01;  /* can do more restrictions */
    }

    if (ispy) {
      nj >>= 1;
      jb = nj-1;
      if (nj & 1)              done |= 0x07;  /* == 3 or 1 */
      else if (nj == 2)        done |= 0x02;  /* can do one more */
    }
    else {
      ja = -((-ja+1)/2) - nu;
      jb = (jb+1)/2 + nu;
      nj = jb - ja + 1;
      if (nj <= omega)         done |= 0x02;  /* can do more restrictions */
    }

    if (ispz) {
      nk >>= 1;
      kb = nk-1;
      if (nk & 1)              done |= 0x07;  /* == 3 or 1 */
      else if (nk == 2)        done |= 0x04;  /* can do one more */
    }
    else {
      ka = -((-ka+1)/2) - nu;
      kb = (kb+1)/2 + nu;
      nk = kb - ka + 1;
      if (nk <= omega)         done |= 0x04;  /* can do more restrictions */
    }

  } while ( ! alldone );
  pm->nlevels = level;

  toplevel = (ispany ? pm->nlevels : pm->nlevels - 1);

  /* ellipsoid axes for grid cutoff weights */
  ni = (int) ceil(2*a/hx) - 1;
  nj = (int) ceil(2*a/hy) - 1;
  nk = (int) ceil(2*a/hz) - 1;
  scaling = 1;  /* initially there is no scaling */
  for (level = 0;  level < toplevel;  level++) {
    if (issprec) {
      p_f = &(pm->gc_f[level]);
      GRID_RESIZE(p_f, float, -ni, 2*ni+1, -nj, 2*nj+1, -nk, 2*nk+1);
    }
    else {
      p = &(pm->gc[level]);
      GRID_RESIZE(p, double, -ni, 2*ni+1, -nj, 2*nj+1, -nk, 2*nk+1);
    }

    if (0 == level) {
      index = 0;
      for (k = -nk;  k <= nk;  k++) {
        for (j = -nj;  j <= nj;  j++) {
          for (i = -ni;  i <= ni;  i++) {
            double s, t, gs, gt, g;
            s = sqrt((i*hx)*(i*hx) + (j*hy)*(j*hy) + (k*hz)*(k*hz)) / a;
            t = 0.5 * s;
            if (t >= 1) {
              g = 0;
            }
            else if (s >= 1) {
              gs = 1/s;
              SPOLY(&gt, &d, t, split);
              g = (gs - 0.5 * gt) / a;
            }
            else {
              SPOLY(&gs, &d, s, split);
              SPOLY(&gt, &d, t, split);
              g = (gs - 0.5 * gt) / a;
            }
            if (issprec) {
              GRID_INDEX_CHECK(p_f, i, j, k);
              ASSERT(p_f->buffer + index == p_f->data + GRID_INDEX(p_f,i,j,k));
              p_f->buffer[index] = (float) g;
            }
            else {
              GRID_INDEX_CHECK(p, i, j, k);
              ASSERT( p->buffer + index == p->data + GRID_INDEX(p, i, j, k) );
              p->buffer[index] = g;
            }
            index++;
          }
        }
      } /* end loops over k-j-i */
    }
    else {
      /* set each level as scaling of h-level */
      if (issprec) {
        const NL_Msmgrid_float *first = &(pm->gc_f[0]);
        index = 0;
        for (k = -nk;  k <= nk;  k++) {
          for (j = -nj;  j <= nj;  j++) {
            for (i = -ni;  i <= ni;  i++) {
              GRID_INDEX_CHECK(p_f, i, j, k);
              ASSERT(p_f->buffer + index == p_f->data + GRID_INDEX(p_f,i,j,k));
              p_f->buffer[index] = (float) (scaling * first->buffer[index]);
              index++;
            }
          }
        } /* for loops */
      } /* if isprec */
      else {
        const NL_Msmgrid_double *first = &(pm->gc[0]);
        index = 0;
        for (k = -nk;  k <= nk;  k++) {
          for (j = -nj;  j <= nj;  j++) {
            for (i = -ni;  i <= ni;  i++) {
              GRID_INDEX_CHECK(p, i, j, k);
              ASSERT( p->buffer + index == p->data + GRID_INDEX(p, i, j, k) );
              p->buffer[index] = scaling * first->buffer[index];
              index++;
            }
          }
        } /* for loops */
      } /* else ! issprec */
    }
    scaling *= 0.5;  /* adjust scaling here to also accommodate toplevel */
  } /* end loop over levels */

  if (toplevel < pm->nlevels) {
    /* nonperiodic in all dimensions,
     * calculate top level weights, ellipsoid axes are length of grid */
    if (issprec) {
      const NL_Msmgrid_float *qhtop_f = &(pm->qh_f[toplevel]);
      ni = qhtop_f->ni - 1;
      nj = qhtop_f->nj - 1;
      nk = qhtop_f->nk - 1;
    }
    else {
      const NL_Msmgrid_double *qhtop = &(pm->qh[toplevel]);
      ni = qhtop->ni - 1;
      nj = qhtop->nj - 1;
      nk = qhtop->nk - 1;
    }
    if (issprec) {
      p_f = &(pm->gc_f[toplevel]);
      GRID_RESIZE(p_f, float, -ni, 2*ni+1, -nj, 2*nj+1, -nk, 2*nk+1);
    }
    else {
      p = &(pm->gc[toplevel]);
      GRID_RESIZE(p, double, -ni, 2*ni+1, -nj, 2*nj+1, -nk, 2*nk+1);
    }
    index = 0;
    for (k = -nk;  k <= nk;  k++) {
      for (j = -nj;  j <= nj;  j++) {
        for (i = -ni;  i <= ni;  i++) {
          double s, gs;
          s = sqrt((i*hx)*(i*hx) + (j*hy)*(j*hy) + (k*hz)*(k*hz)) / a;
          if (s >= 1) {
            gs = 1/s;
          }
          else {
            SPOLY(&gs, &d, s, split);
          }
          if (issprec) {
            GRID_INDEX_CHECK(p_f, i, j, k);
            ASSERT(p_f->buffer + index == p_f->data + GRID_INDEX(p_f,i,j,k));
            p_f->buffer[index] = (float) (scaling * gs/a);
          }
          else {
            GRID_INDEX_CHECK(p, i, j, k);
            ASSERT( p->buffer + index == p->data + GRID_INDEX(p, i, j, k) );
            p->buffer[index] = scaling * gs/a;
          }
          index++;
        }
      }
    } /* end loops over k-j-i for coarsest level weights */
  }

  /* calculate self energy factor for splitting */
  if (1) {
    double s, gs;
    s = 0;
    SPOLY(&gs, &d, s, split);
    if (issprec) {
      pm->gzero_f = (float) (gs/a);
    }
    else {
      pm->gzero = gs/a;
    }
  }

  return NL_MSM_SUCCESS;
}
