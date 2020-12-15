#include <stdlib.h>
#include "sYlm.h"



// pi
const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068;


static inline double fac(const int n)
{
   double result = 1;
   
   for (int i=2; i <= n; i++)
      result *= i;
	    
   return result;
}


// coefficient function
static inline double Cslm(const int s, const int l, const int m)
{
   assert(l >= 0);

   return sqrt( l*l * (4.0*l*l - 1.0) / ( (l*l - m*m) * (l*l - s*s) ) );
}



// recursion function
double s_lambda_lm(const int s, const int l, const int m, const double x)
{
   double Pm;

   assert(m >= 0);
   assert(l >= 0);
   assert(m >= s);

   Pm = pow(-0.5, m);
/*
   if (l == m) // recursion end
      return pow(-2.0, -m) * sqrt( fac(2*m + 1) * 1.0 / ( 4.0*pi * fac(m+s) * fac(m-s) ) * pow(1.0-x, m+s) * pow(1.0+x, m-s) );
   else if (l >= m+2)
      return  (x + s*m * 1.0 / ( l * (l-1.0) ) ) * Cslm(s, l, m) * s_lambda_lm(s, l-1, m, x) 
             - Cslm(s, l, m) * 1.0 / Cslm(s, l-1, m) * s_lambda_lm(s, l-2, m, x);
   else if (l >= m+1)
      return  (x + s*m * 1.0 / ( l * (l-1.0) ) ) * Cslm(s, l, m) * s_lambda_lm(s, l-1, m, x);
   else
      return 0.0;
*/

   if (m !=  s) Pm = Pm * pow(1.0+x, (m-s)*1.0/2);
   if (m != -s) Pm = Pm * pow(1.0-x, (m+s)*1.0/2);
   
   Pm = Pm * sqrt( fac(2*m + 1) * 1.0 / ( 4.0*pi * fac(m+s) * fac(m-s) ) );
   
   if (l == m)
      return Pm;
   
   double Pm1 = (x + s*1.0/(m+1) ) * Cslm(s, m+1, m) * Pm;
   
   if (l == m+1)
      return Pm1;
   else
   {
      double Pn;
      for (int n=m+2; n <= l; n++)
      {
         Pn = (x + s*m * 1.0 / ( n * (n-1.0) ) ) * Cslm(s, n, m) * Pm1 
             - Cslm(s, n, m) * 1.0 / Cslm(s, n-1, m) * Pm;
         Pm = Pm1;
         Pm1 = Pn;
         
      }
      return Pn;
   }
}





void sYlm(const int ss, const int ll, const int mm, const double theta, const double phi, double* const res_re, double* const res_im)
{
   int s, l, m;
   double Pm = 1.0;

   l = ll;
   m = mm;
   s = ss;

   *res_re = 0.0;
   *res_im = 0.0;

   if (l < 0)
      return;
   if (abs(m) > l || l < abs(s))
      return;

   if (abs(mm) < abs(ss))
   {
      s=mm;
      m=ss;
      if ((m+s) % 2)
         Pm  = -Pm;
   }
   
   if (m < 0)
   {
      s=-s;
      m=-m;
      if ((m+s) % 2)
         Pm  = -Pm;
   }

   double result = Pm * s_lambda_lm(s, l, m, cos(theta));

   if (!finite(result)) result = 0.0;

   *res_re = result * cos(mm*phi);
   *res_im = result * sin(mm*phi);
}



