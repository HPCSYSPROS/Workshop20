fp t35;
fp t41;
fp t38;
fp t28;
fp t37;
fp t40;
fp t36;
fp t39;
fp t34;
fp t33;
fp t30;
fp t29;
fp t27;
      t35 = RATIONAL(1.0,2.0);
      t41 = x*t35;
      t38 = y*y;
      t28 = t38*t41;
      t37 = RATIONAL(-1.0,4.0);
      t40 = t28+t37*y;
      t36 = RATIONAL(1.0,4.0);
      t39 = t28+t36*y;
      t34 = RATIONAL(-1.0,2.0);
      t33 = t37*t38;
      t30 = t36*t38;
      t29 = t34*x*y;
      t27 = y*t41;
      coeffs_dx->coeff_m1_m1 = t29+t33+t39;
      coeffs_dx->coeff_0_m1 = (-t38+y)*x;
      coeffs_dx->coeff_p1_m1 = t29+t30+t40;
      coeffs_dx->coeff_m1_0 = t34+x+(t35-x)*t38;
      coeffs_dx->coeff_0_0 = (RATIONAL(-2.0,1.0)+RATIONAL(2.0,1.0)*t38)*x;
      coeffs_dx->coeff_p1_0 = t35+x+(-x+t34)*t38;
      coeffs_dx->coeff_m1_p1 = t27+t33+t40;
      coeffs_dx->coeff_0_p1 = (-y-t38)*x;
      coeffs_dx->coeff_p1_p1 = t27+t30+t39;
