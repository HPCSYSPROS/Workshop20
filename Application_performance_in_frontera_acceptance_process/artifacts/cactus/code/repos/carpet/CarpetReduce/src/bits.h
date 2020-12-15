#ifndef BITS_H
#define BITS_H

/* a mask with exactly bit n set */
inline unsigned BMSK(unsigned n) { return 1U << n; }

/* whether bit n in value v is set */
inline unsigned BGET(unsigned v, unsigned n) { return !!(v & BMSK(n)); }

/* set, clear, or invert bit n in value v */
inline unsigned BSET(unsigned v, unsigned n) { return v | BMSK(n); }
inline unsigned BCLR(unsigned v, unsigned n) { return v & ~BMSK(n); }
inline unsigned BINV(unsigned v, unsigned n) { return v ^ BMSK(n); }

/* set bit n in value v to b */
inline unsigned BCPY(unsigned v, unsigned n, unsigned b) {
  return BCLR(v, n) | (-!!b & BMSK(n));
}

/* count the number of set bits */
inline unsigned BCNT(unsigned v) {
  v = (v & 0x55555555U) + ((v >> 1) & 0x55555555U);
  v = (v & 0x33333333U) + ((v >> 2) & 0x33333333U);
  v = (v + (v >> 4)) & 0x0f0f0f0fU;
  v = v + (v >> 8);
  v = v + (v >> 16);
  return v & 0xff;
}

#endif /* #ifndef BITS_H */
