/* Add the contributions to T_munu that we gathered via the AddToTmunu
   mechanism, using the CalcTmunu.inc mechanism.  */



/* Only add if we are not included from this thorn; otherwise we would
   be adding to ourselves.  */
#ifndef TMUNUBASE_SETTMUNU

#ifdef CCODE

  if (* stress_energy_2_state)
  {

#ifndef NDEBUG
    if ((ijk) < 0 || (ijk) >= cctk_ash[0] * cctk_ash[1] * cctk_ash[2])
    {
      CCTK_WARN (0, "Illegal array index");
    }
#endif

    Ttt += eT2tt[ijk];

    Ttx += eT2tx[ijk];
    Tty += eT2ty[ijk];
    Ttz += eT2tz[ijk];

    Txx += eT2xx[ijk];
    Txy += eT2xy[ijk];
    Txz += eT2xz[ijk];
    Tyy += eT2yy[ijk];
    Tyz += eT2yz[ijk];
    Tzz += eT2zz[ijk];

  }

#endif

#ifdef FCODE

      if (stress_energy_2_state .ne. 0) then

        Ttt = Ttt + eT2tt(i,j,k)
                 
        Ttx = Ttx + eT2tx(i,j,k)
        Tty = Tty + eT2ty(i,j,k)
        Ttz = Ttz + eT2tz(i,j,k)
                 
        Txx = Txx + eT2xx(i,j,k)
        Txy = Txy + eT2xy(i,j,k)
        Txz = Txz + eT2xz(i,j,k)
        Tyy = Tyy + eT2yy(i,j,k)
        Tyz = Tyz + eT2yz(i,j,k)
        Tzz = Tzz + eT2zz(i,j,k)

      end if

#endif

#endif /* ifndef TMUNUBASE_SETTMUNU */
