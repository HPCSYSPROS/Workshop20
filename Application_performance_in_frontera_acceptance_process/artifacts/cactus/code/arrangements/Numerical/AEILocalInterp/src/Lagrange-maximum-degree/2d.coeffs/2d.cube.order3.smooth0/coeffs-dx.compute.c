fp t199;
fp t185;
fp t178;
fp t181;
fp t186;
fp t173;
fp t198;
fp t159;
fp t197;
fp t182;
fp t175;
fp t196;
fp t161;
fp t195;
fp t162;
fp t194;
fp t160;
fp t193;
fp t179;
fp t184;
fp t176;
fp t192;
fp t183;
fp t171;
fp t191;
fp t180;
fp t190;
fp t189;
fp t177;
fp t188;
fp t187;
fp t174;
fp t169;
fp t167;
fp t165;
fp t164;
fp t163;
      t199 = x*y;
      t185 = x*x;
      t178 = RATIONAL(-3.0,8.0)*t185;
      t181 = RATIONAL(1.0,40.0);
      t186 = y*y;
      t173 = t181*t186;
      t198 = t178+t173;
      t159 = RATIONAL(-3.0,20.0)*t199;
      t197 = RATIONAL(13.0,40.0)*x+t159;
      t182 = RATIONAL(-1.0,40.0);
      t175 = t182*t186;
      t196 = t178+t175;
      t161 = RATIONAL(-1.0,20.0)*t199;
      t195 = t161+RATIONAL(11.0,40.0)*x;
      t162 = RATIONAL(3.0,20.0)*t199;
      t194 = RATIONAL(7.0,40.0)*x+t162;
      t160 = RATIONAL(1.0,20.0)*t199;
      t193 = RATIONAL(9.0,40.0)*x+t160;
      t179 = RATIONAL(-1.0,8.0)*t185;
      t184 = RATIONAL(-3.0,40.0);
      t176 = t184*t186;
      t192 = t179+t176;
      t183 = RATIONAL(3.0,40.0);
      t171 = t183*t186;
      t191 = t179+t171;
      t180 = RATIONAL(3.0,8.0)*t185;
      t190 = t180+t175;
      t189 = t180+t173;
      t177 = RATIONAL(1.0,8.0)*t185;
      t188 = t171+t177;
      t187 = t176+t177;
      t174 = RATIONAL(-1.0,50.0)*y;
      t169 = RATIONAL(2.0,25.0)*y;
      t167 = RATIONAL(-9.0,100.0)*y;
      t165 = RATIONAL(-1.0,100.0)*y;
      t164 = RATIONAL(7.0,100.0)*y;
      t163 = RATIONAL(-13.0,100.0)*y;
      coeffs_dx->coeff_m1_m1 = RATIONAL(6.0,25.0)*y+RATIONAL(-109.0,1200.0)+
t192+t197;
      coeffs_dx->coeff_0_m1 = t162+t174+RATIONAL(-31.0,400.0)+RATIONAL(-23.0,
40.0)*x+t190;
      coeffs_dx->coeff_p1_m1 = RATIONAL(111.0,400.0)+t163+t194+t198;
      coeffs_dx->coeff_p2_m1 = t183*x+RATIONAL(-131.0,1200.0)+t167+t159+t188;
      coeffs_dx->coeff_m1_0 = RATIONAL(-223.0,1200.0)+t174+t191+t195;
      coeffs_dx->coeff_0_0 = RATIONAL(-1.0,25.0)*y+t160+RATIONAL(-21.0,40.0)*x+
RATIONAL(-57.0,400.0)+t189;
      coeffs_dx->coeff_p1_0 = t165+RATIONAL(117.0,400.0)+t193+t196;
      coeffs_dx->coeff_p2_0 = t161+t181*x+RATIONAL(43.0,1200.0)+t164+t187;
      coeffs_dx->coeff_m1_p1 = RATIONAL(-157.0,1200.0)+t163+t191+t193;
      coeffs_dx->coeff_0_p1 = RATIONAL(-63.0,400.0)+t165+t161+RATIONAL(-19.0,
40.0)*x+t189;
      coeffs_dx->coeff_p1_p1 = RATIONAL(103.0,400.0)+RATIONAL(3.0,50.0)*y+t195+
t196;
      coeffs_dx->coeff_p2_p1 = t160+t169+RATIONAL(37.0,1200.0)+t182*x+t187;
      coeffs_dx->coeff_m1_p2 = RATIONAL(89.0,1200.0)+t167+t192+t194;
      coeffs_dx->coeff_0_p2 = t159+t164+RATIONAL(-17.0,40.0)*x+RATIONAL(-49.0,
400.0)+t190;
      coeffs_dx->coeff_p1_p2 = RATIONAL(69.0,400.0)+t169+t197+t198;
      coeffs_dx->coeff_p2_p2 = RATIONAL(-149.0,1200.0)+t162+RATIONAL(-3.0,50.0)
*y+t184*x+t188;
