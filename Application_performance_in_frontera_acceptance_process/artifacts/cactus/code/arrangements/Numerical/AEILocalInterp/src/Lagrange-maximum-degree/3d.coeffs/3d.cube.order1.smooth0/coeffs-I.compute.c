fp t8;
fp t2;
fp t7;
fp t4;
fp t13;
fp t1;
fp t6;
fp t12;
fp t11;
fp t10;
fp t5;
fp t9;
fp t3;
      t8 = RATIONAL(1.0,4.0);
      t2 = t8*y;
      t7 = RATIONAL(-1.0,4.0);
      t4 = t7*z;
      t13 = t2+t4;
      t1 = t8*z;
      t6 = t7*y;
      t12 = t1+t6;
      t11 = t1+t2;
      t10 = t4+t6;
      t5 = t7*x;
      t9 = t5+t8;
      t3 = t8*x;
      coeffs_I->coeff_0_0_0 = t5+RATIONAL(1.0,2.0)+t10;
      coeffs_I->coeff_p1_0_0 = t8+t3+t10;
      coeffs_I->coeff_0_p1_0 = t9+t13;
      coeffs_I->coeff_p1_p1_0 = t3+t13;
      coeffs_I->coeff_0_0_p1 = t9+t12;
      coeffs_I->coeff_p1_0_p1 = t3+t12;
      coeffs_I->coeff_0_p1_p1 = t5+t11;
      coeffs_I->coeff_p1_p1_p1 = t3+t7+t11;
