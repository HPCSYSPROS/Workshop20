fp t704;
fp t708;
fp t685;
fp t709;
fp t688;
fp t726;
fp t683;
fp t695;
fp t725;
fp t724;
fp t723;
fp t722;
fp t679;
fp t692;
fp t694;
fp t687;
fp t681;
fp t721;
fp t693;
fp t720;
fp t686;
fp t719;
fp t718;
fp t682;
fp t717;
fp t678;
fp t716;
fp t684;
fp t715;
fp t714;
fp t713;
fp t712;
fp t689;
fp t711;
fp t710;
fp t677;
fp t676;
fp t675;
fp t674;
fp t673;
fp t672;
      t704 = RATIONAL(1.0,10.0);
      t708 = x*x;
      t685 = t704*t708;
      t709 = y*y;
      t688 = RATIONAL(2.0,49.0)*t709;
      t726 = t685+t688+RATIONAL(-289.0,2940.0);
      t683 = RATIONAL(-1.0,49.0)*t709;
      t695 = RATIONAL(-2.0,5.0)*t708;
      t725 = t683+t695+RATIONAL(226.0,735.0);
      t724 = RATIONAL(1.0,98.0)*t709+t695+RATIONAL(181.0,735.0);
      t723 = t683+t685+RATIONAL(71.0,2940.0);
      t722 = x*y;
      t679 = t704*x;
      t692 = RATIONAL(-1.0,5.0)*x;
      t694 = RATIONAL(1.0,5.0)*x;
      t687 = RATIONAL(-1.0,10.0)*x;
      t681 = RATIONAL(-2.0,49.0)*t709;
      t721 = t681+t685+RATIONAL(191.0,2940.0);
      t693 = RATIONAL(3.0,5.0)*t708;
      t720 = t681+t693+RATIONAL(-41.0,98.0);
      t686 = RATIONAL(1.0,49.0)*t709;
      t719 = t686+t695+RATIONAL(166.0,735.0);
      t718 = t686+t693+RATIONAL(-53.0,98.0);
      t682 = RATIONAL(2.0,35.0)*y;
      t717 = t682+t726;
      t678 = RATIONAL(-1.0,35.0)*y;
      t716 = t678+t725;
      t684 = RATIONAL(1.0,35.0)*y;
      t715 = t684+t725;
      t714 = t679+t723;
      t713 = t687+t723;
      t712 = RATIONAL(-1.0,70.0)*y+t724;
      t689 = RATIONAL(-2.0,35.0)*y;
      t711 = t689+t726;
      t710 = RATIONAL(1.0,70.0)*y+t724;
      t677 = y*t692;
      t676 = y*t694;
      t675 = RATIONAL(1.0,20.0)*t722;
      t674 = y*t687;
      t673 = RATIONAL(-1.0,20.0)*t722;
      t672 = y*t679;
      coeffs_dxx->coeff_m2_m2 = t672+t687+t711;
      coeffs_dxx->coeff_m1_m2 = t677+t694+t715;
      coeffs_dxx->coeff_0_m2 = t682+t720;
      coeffs_dxx->coeff_p1_m2 = t676+t692+t715;
      coeffs_dxx->coeff_p2_m2 = t674+t679+t711;
      coeffs_dxx->coeff_m2_m1 = t675+t678+t713;
      coeffs_dxx->coeff_m1_m1 = t694+t674+t710;
      coeffs_dxx->coeff_0_m1 = t684+t718;
      coeffs_dxx->coeff_p1_m1 = t692+t672+t710;
      coeffs_dxx->coeff_p2_m1 = t678+t673+t714;
      coeffs_dxx->coeff_m2_0 = t687+t721;
      coeffs_dxx->coeff_m1_0 = t694+t719;
      coeffs_dxx->coeff_0_0 = t693+RATIONAL(-57.0,98.0)+t688;
      coeffs_dxx->coeff_p1_0 = t692+t719;
      coeffs_dxx->coeff_p2_0 = t679+t721;
      coeffs_dxx->coeff_m2_p1 = t684+t673+t713;
      coeffs_dxx->coeff_m1_p1 = t672+t694+t712;
      coeffs_dxx->coeff_0_p1 = t678+t718;
      coeffs_dxx->coeff_p1_p1 = t692+t674+t712;
      coeffs_dxx->coeff_p2_p1 = t684+t675+t714;
      coeffs_dxx->coeff_m2_p2 = t687+t674+t717;
      coeffs_dxx->coeff_m1_p2 = t694+t676+t716;
      coeffs_dxx->coeff_0_p2 = t689+t720;
      coeffs_dxx->coeff_p1_p2 = t677+t692+t716;
      coeffs_dxx->coeff_p2_p2 = t679+t672+t717;