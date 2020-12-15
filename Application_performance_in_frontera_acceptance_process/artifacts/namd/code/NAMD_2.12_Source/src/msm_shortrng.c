/** msm_shortrng.c
 *
 * Compute the short-range part of MSM forces.
 *
 * Perform spatial hashing of atoms into bins.  The implementation
 * uses the simplest approach where the bin lengths are at least as
 * large as the cutoff distance, so only nearest neighbor bins need
 * be checked.  The bins are represented by a simple linked list of
 * atom indices using fixed array storage: an array of int for the
 * bins giving the index of a "head" atom in the bin and an array of
 * int for the "next" atom in the bin, where -1 denotes the end of
 * the list.
 *
 * The formulation of bins is determined for an arbitrary set of
 * basis vectors and the spatial hashing is performed in reciprocal
 * space, so this part of the MSM code is ready to support
 * non-orthogonal basis vectors.
 *
 * XXX The virial is not calculated.
 *
 * XXX The idea of excluded interaction pairs is totally disregarded.
 * This should be fine for calculating normal systems of atoms, where
 * the full 1/r interaction is later subtracted out for excluded pairs
 * without causing excessive numerical error.  However, doing so will
 * be a problem for systems containing Drude particles, possibly also
 * for lone pairs, and is unsuitable for calculating Lennard-Jones.
 */

#include "msm_defn.h"

static int setup_bin_data(NL_Msm *pm);
static int spatial_hashing(NL_Msm *pm);
static int bin_evaluation_1away(NL_Msm *pm);
static int bin_evaluation_k_away(NL_Msm *pm);


int NL_msm_compute_short_range(NL_Msm *pm) {
  int rc = 0;  /* return code */

  if (pm->maxatoms < pm->numatoms) {
    rc = setup_bin_data(pm);
    if (rc) return rc;
  }
  rc = spatial_hashing(pm);
  if (rc) return rc;
  if (pm->msmflags & NL_MSM_COMPUTE_1AWAY) {
    rc = bin_evaluation_1away(pm);
  }
  else {
    rc = bin_evaluation_k_away(pm);
  }
  if (rc) return rc;
  return NL_MSM_SUCCESS;
}


int setup_bin_data(NL_Msm *pm) {
  void *vptr = NULL;

  ASSERT(pm->maxatoms < pm->numatoms);
  vptr = malloc(pm->numatoms * sizeof(int));
  if (vptr == NULL) return NL_MSM_ERROR_MALLOC;
  pm->next = (int *) vptr;
  pm->maxatoms = pm->numatoms;
  return NL_MSM_SUCCESS;
}


int spatial_hashing(NL_Msm *pm) {
  double d[3], s[3];
  const double *atom = pm->atom;  /* stored x/y/z/q */
  const double *c = pm->cellcenter;
  const double *ru = pm->recipvec1;
  const double *rv = pm->recipvec2;
  const double *rw = pm->recipvec3;
  int *bin = pm->bin;
  int *next = pm->next;
  int numbins = pm->numbins;
  int nbx = pm->nbx;
  int nby = pm->nby;
  int nbz = pm->nbz;
  int numatoms = pm->numatoms;
  int n, ib, jb, kb, index;

  /* reset all bin and next IDs */
  for (n = 0;  n < numbins;  n++)  bin[n] = -1;
  for (n = 0;  n < numatoms;  n++)  next[n] = -1;

  for (n = 0;  n < numatoms;  n++, atom += 4) {
    /* transform coordinates to reciprocal space */
    d[X] = atom[X] - c[X];
    d[Y] = atom[Y] - c[Y];
    d[Z] = atom[Z] - c[Z];
    s[X] = ru[X] * d[X] + ru[Y] * d[Y] + ru[Z] * d[Z];
    s[Y] = rv[X] * d[X] + rv[Y] * d[Y] + rv[Z] * d[Z];
    s[Z] = rw[X] * d[X] + rw[Y] * d[Y] + rw[Z] * d[Z];
    /* determine bin indexing in 3D */
    ib = (int) floor((s[X] + 0.5) * nbx);
    jb = (int) floor((s[Y] + 0.5) * nby);
    kb = (int) floor((s[Z] + 0.5) * nbz);
    /* assume coordinate is within defined cell;
     * adjust bin index in case of roundoff error */
    if      (ib < 0)    ib = 0;
    else if (ib >= nbx) ib = nbx - 1;
    if      (jb < 0)    jb = 0;
    else if (jb >= nby) jb = nby - 1;
    if      (kb < 0)    kb = 0;
    else if (kb >= nbz) kb = nbz - 1;
    index = (kb*nby + jb)*nbx + ib;  /* flatten 3D indexing into 1D */
    next[n] = bin[index];  /* attach atom index to front of linked list */
    bin[index] = n;
  }
  return NL_MSM_SUCCESS;
}


int bin_evaluation_1away(NL_Msm *pm) {
  enum { NUM_NBRBINS = 14 };
  int nbrbin[3*NUM_NBRBINS] = {
    0,0,0,  1,0,0,  -1,1,0,  0,1,0,  1,1,0,  -1,-1,1,  0,-1,1,
    1,-1,1,  -1,0,1,  0,0,1,  1,0,1,  -1,1,1,  0,1,1,  1,1,1
  };
  int nbx = pm->nbx;
  int nby = pm->nby;
  int nbz = pm->nbz;
  int aindex, bindex;
  int ia, ja, ka, n, i, j;
  int ispx = ((pm->msmflags & NL_MSM_PERIODIC_VEC1) != 0);
  int ispy = ((pm->msmflags & NL_MSM_PERIODIC_VEC2) != 0);
  int ispz = ((pm->msmflags & NL_MSM_PERIODIC_VEC3) != 0);
  int split = pm->split;
  const double *u = pm->cellvec1;
  const double *v = pm->cellvec2;
  const double *w = pm->cellvec3;
  const int *bin = pm->bin;
  const int *next = pm->next;
  const double *atom = pm->atom;
  double *force = pm->felec;
  double a2 = pm->a * pm->a;  /* cutoff^2 */
  double a_1 = 1 / pm->a;     /* 1 / cutoff */
  double a_2 = a_1 * a_1;     /* 1 / cutoff^2 */
  double u_elec = 0;  /* accumulate potential energy from electrostatics */

  for (aindex = 0, ka = 0;  ka < nbz;  ka++) {
    for (ja = 0;  ja < nby;  ja++) {
      for (ia = 0;  ia < nbx;  ia++, aindex++) {  /* loop over bins A */

        for (n = 0;  n < NUM_NBRBINS;  n++) { /* loop B-bin neighborhood of A */

          double p[3] = { 0,0,0 };  /* periodic wrapping vector for B bin */

          int ib = ia + nbrbin[3*n + X];  /* index for B bin */
          int jb = ja + nbrbin[3*n + Y];  /* index for B bin */
          int kb = ka + nbrbin[3*n + Z];  /* index for B bin */

          /* do wrap around for bin index outside of periodic dimension range,
           * short-circuit loop for bin index outside non-periodic range */
          if (ispx) {
            if (ib < 0) {
              ib += nbx;  p[X] -= u[X];  p[Y] -= u[Y];  p[Z] -= u[Z];
            }
            else if (ib >= nbx) {
              ib -= nbx;  p[X] += u[X];  p[Y] += u[Y];  p[Z] += u[Z];
            }
          }
          else if (ib < 0 || ib >= nbx) continue;

          if (ispy) {
            if (jb < 0) {
              jb += nby;  p[X] -= v[X];  p[Y] -= v[Y];  p[Z] -= v[Z];
            }
            else if (jb >= nby) {
              jb -= nby;  p[X] += v[X];  p[Y] += v[Y];  p[Z] += v[Z];
            }
          }
          else if (jb < 0 || jb >= nby) continue;

          if (ispz) {
            if (kb < 0) {
              kb += nbz;  p[X] -= w[X];  p[Y] -= w[Y];  p[Z] -= w[Z];
            }
            else if (kb >= nbz) {
              kb -= nbz;  p[X] += w[X];  p[Y] += w[Y];  p[Z] += w[Z];
            }
          }
          else if (kb < 0 || kb >= nbz) continue;

          /* flat 1D index for B bin, after doing wrap around */
          bindex = (kb*nby + jb)*nbx + ib;

          for (j = bin[bindex];  j != -1;  j = next[j]) { /* loop over B bin */
            double force_j[3] = { 0,0,0 };
            double atom_j[3];
            double qj = atom[4*j + Q];
            int ihead;

            atom_j[X] = atom[4*j + X] + p[X];
            atom_j[Y] = atom[4*j + Y] + p[Y];
            atom_j[Z] = atom[4*j + Z] + p[Z];

            ihead = bin[aindex];
            /* for self bin (A==B) interactions, visit each pair only once */
            if (n == 0) ihead = next[j];  /* i.e. aindex == bindex */

            for (i = ihead;  i != -1;  i = next[i]) { /* loop over A bin */
              double rij[3];
              double r2;

              rij[X] = atom_j[X] - atom[4*i + X];
              rij[Y] = atom_j[Y] - atom[4*i + Y];
              rij[Z] = atom_j[Z] - atom[4*i + Z];

              r2 = rij[X] * rij[X] + rij[Y] * rij[Y] + rij[Z] * rij[Z];

              if (r2 < a2) {
                double fij[3];
                double qq = qj * atom[4*i + Q];  /* combined charge */
                double r;      /* length of vector r_ij */
                double r_1;    /* 1/r */
                double r_2;    /* 1/r^2 */
                double r_a;    /* r/a */
                double g;      /* normalized smoothing g(R), R=r/a */
                double dg;     /* (d/dR)g(R) */
                double ue;     /* U_elec(r) */
                double due_r;  /* (1/r)*(d/dr)U_elec(r) */

                r = sqrt(r2);
                r_1 = 1/r;
                r_2 = r_1 * r_1;

                /* calculate MSM splitting */
                r_a = r * a_1;
                SPOLY(&g, &dg, r_a, split);

                ue = qq * (r_1 - a_1 * g);
                due_r = qq * r_1 * (-r_2 - a_2 * dg);

                fij[X] = -rij[X] * due_r;
                fij[Y] = -rij[Y] * due_r;
                fij[Z] = -rij[Z] * due_r;
                force[3*i + X] -= fij[X];
                force[3*i + Y] -= fij[Y];
                force[3*i + Z] -= fij[Z];
                force_j[X] += fij[X];
                force_j[Y] += fij[Y];
                force_j[Z] += fij[Z];

                u_elec += ue;

                /* XXX virial? */

              } /* end if r2 < cutoff2 */

            } /* end loop over A bin */

            force[3*j + X] += force_j[X];
            force[3*j + Y] += force_j[Y];
            force[3*j + Z] += force_j[Z];
          } /* end loop over B bin */

        } /* end loop B-bin neighborhood of A */

      } /* end loop over bins A */
    }
  }

  pm->uelec += u_elec;

  return NL_MSM_SUCCESS;
}


int bin_evaluation_k_away(NL_Msm *pm) {
  int *nbrhd = pm->nbrhd;
  int nbrhdlen = pm->nbrhdlen;
  int nbx = pm->nbx;
  int nby = pm->nby;
  int nbz = pm->nbz;
  int aindex, bindex;
  int ia, ja, ka, n, i, j;
  int ispx = ((pm->msmflags & NL_MSM_PERIODIC_VEC1) != 0);
  int ispy = ((pm->msmflags & NL_MSM_PERIODIC_VEC2) != 0);
  int ispz = ((pm->msmflags & NL_MSM_PERIODIC_VEC3) != 0);
  int split = pm->split;
  const double *u = pm->cellvec1;
  const double *v = pm->cellvec2;
  const double *w = pm->cellvec3;
  const int *bin = pm->bin;
  const int *next = pm->next;
  const double *atom = pm->atom;
  double *force = pm->felec;
  double a2 = pm->a * pm->a;  /* cutoff^2 */
  double a_1 = 1 / pm->a;     /* 1 / cutoff */
  double a_2 = a_1 * a_1;     /* 1 / cutoff^2 */
  double u_elec = 0;  /* accumulate potential energy from electrostatics */

  for (aindex = 0, ka = 0;  ka < nbz;  ka++) {
    for (ja = 0;  ja < nby;  ja++) {
      for (ia = 0;  ia < nbx;  ia++, aindex++) {  /* loop over bins A */

        for (n = 0;  n < nbrhdlen;  n++) { /* loop B-bin neighborhood of A */

          double p[3] = { 0,0,0 };  /* periodic wrapping vector for B bin */

          int ib = ia + nbrhd[3*n + X];  /* index for B bin */
          int jb = ja + nbrhd[3*n + Y];  /* index for B bin */
          int kb = ka + nbrhd[3*n + Z];  /* index for B bin */

          /* do wrap around for bin index outside of periodic dimension range,
           * short-circuit loop for bin index outside non-periodic range */
          if (ispx) {
            if (ib < 0) {
              ib += nbx;  p[X] -= u[X];  p[Y] -= u[Y];  p[Z] -= u[Z];
            }
            else if (ib >= nbx) {
              ib -= nbx;  p[X] += u[X];  p[Y] += u[Y];  p[Z] += u[Z];
            }
          }
          else if (ib < 0 || ib >= nbx) continue;

          if (ispy) {
            if (jb < 0) {
              jb += nby;  p[X] -= v[X];  p[Y] -= v[Y];  p[Z] -= v[Z];
            }
            else if (jb >= nby) {
              jb -= nby;  p[X] += v[X];  p[Y] += v[Y];  p[Z] += v[Z];
            }
          }
          else if (jb < 0 || jb >= nby) continue;

          if (ispz) {
            if (kb < 0) {
              kb += nbz;  p[X] -= w[X];  p[Y] -= w[Y];  p[Z] -= w[Z];
            }
            else if (kb >= nbz) {
              kb -= nbz;  p[X] += w[X];  p[Y] += w[Y];  p[Z] += w[Z];
            }
          }
          else if (kb < 0 || kb >= nbz) continue;

          /* flat 1D index for B bin, after doing wrap around */
          bindex = (kb*nby + jb)*nbx + ib;

          for (j = bin[bindex];  j != -1;  j = next[j]) { /* loop over B bin */
            double force_j[3] = { 0,0,0 };
            double atom_j[3];
            double qj = atom[4*j + Q];
            int ihead;

            atom_j[X] = atom[4*j + X] + p[X];
            atom_j[Y] = atom[4*j + Y] + p[Y];
            atom_j[Z] = atom[4*j + Z] + p[Z];

            ihead = bin[aindex];
            /* for self bin (A==B) interactions, visit each pair only once */
            if (n == 0) ihead = next[j];  /* i.e. aindex == bindex */

            for (i = ihead;  i != -1;  i = next[i]) { /* loop over A bin */
              double rij[3];
              double r2;

              rij[X] = atom_j[X] - atom[4*i + X];
              rij[Y] = atom_j[Y] - atom[4*i + Y];
              rij[Z] = atom_j[Z] - atom[4*i + Z];

              r2 = rij[X] * rij[X] + rij[Y] * rij[Y] + rij[Z] * rij[Z];

              if (r2 < a2) {
                double fij[3];
                double qq = qj * atom[4*i + Q];  /* combined charge */
                double r;      /* length of vector r_ij */
                double r_1;    /* 1/r */
                double r_2;    /* 1/r^2 */
                double r_a;    /* r/a */
                double g;      /* normalized smoothing g(R), R=r/a */
                double dg;     /* (d/dR)g(R) */
                double ue;     /* U_elec(r) */
                double due_r;  /* (1/r)*(d/dr)U_elec(r) */

                r = sqrt(r2);
                r_1 = 1/r;
                r_2 = r_1 * r_1;

                /* calculate MSM splitting */
                r_a = r * a_1;
                SPOLY(&g, &dg, r_a, split);

                ue = qq * (r_1 - a_1 * g);
                due_r = qq * r_1 * (-r_2 - a_2 * dg);

                fij[X] = -rij[X] * due_r;
                fij[Y] = -rij[Y] * due_r;
                fij[Z] = -rij[Z] * due_r;
                force[3*i + X] -= fij[X];
                force[3*i + Y] -= fij[Y];
                force[3*i + Z] -= fij[Z];
                force_j[X] += fij[X];
                force_j[Y] += fij[Y];
                force_j[Z] += fij[Z];

                u_elec += ue;

                /* XXX virial? */

              } /* end if r2 < cutoff2 */

            } /* end loop over A bin */

            force[3*j + X] += force_j[X];
            force[3*j + Y] += force_j[Y];
            force[3*j + Z] += force_j[Z];
          } /* end loop over B bin */

        } /* end loop B-bin neighborhood of A */

      } /* end loop over bins A */
    }
  }

  pm->uelec += u_elec;

  return NL_MSM_SUCCESS;
}
