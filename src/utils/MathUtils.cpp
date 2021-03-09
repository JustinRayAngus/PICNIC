#include "MathUtils.H"
//#include "MathUtilsF_F.H"

#include <cmath>
#include <cstdlib>

#include "NamespaceHeader.H"
      
double MathUtils::errorfun( const double&  a_x )
{
   // errorfun(x) = 2/sqrt(pi)*int_0^x(exp(-x^2))dx
   //
   double soln = erf(a_x);
   return soln;

}

double MathUtils::errorinv( const double&  a_x )
{
   CH_assert( abs(a_x)<=1.0 );

   double soln;
   double z;

   double y0 = 0.7;
   vector<double> a{0.886226899,  -1.645349621,  0.914624893, -0.140543331}; 
   vector<double> b{-2.118377725,  1.442710462, -0.329097515,  0.012229801};
   vector<double> c{-1.970840454, -1.624906493,  3.429567803,  1.641345311}; 
   vector<double> d{ 3.543889200,  1.637067800};

   if(abs(a_x)==1.0) {
      soln = -a_x*log(0.0);
   }
   else if(a_x<-y0) {
      z = sqrt(-log((1.0+a_x)/2.0));
      soln = -(((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1.0);  
   } 
   else {
      if(a_x<y0) {
         z = a_x*a_x;
         soln = a_x*(((a[3]*z+a[2])*z+a[1])*z+a[0])/((((b[3]*z+b[2])*z+b[1])*z+b[0])*z+1.0);
      }
      else {
         z = sqrt(-log((1.0-a_x)/2.0));
         soln = (((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1.0);
      }
      soln = soln - (erf(soln) - a_x)/(2.0/sqrt(M_PI)*exp(-soln*soln));
      soln = soln - (erf(soln) - a_x)/(2.0/sqrt(M_PI)*exp(-soln*soln));
   }

   return soln;

}

double MathUtils::rand()
{
   
   static std::random_device rd;
   static std::mt19937 gen(rd());
   static std::uniform_real_distribution<> dis(0,1);
   double randNum = dis(gen);

   return randNum;

}

int MathUtils::randInt(const int& A, const int& B)
{
   
   static std::random_device randDev;
   static std::mt19937 twister(randDev());
   static std::uniform_int_distribution<int> dist;

   dist.param(std::uniform_int_distribution<int>::param_type(A,B));

   return dist(twister);

}





