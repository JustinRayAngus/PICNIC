#include "MathUtils.H"
//#include "MathUtilsF_F.H"

#include <cmath>
#include <cstdlib>

#include "NamespaceHeader.H"

void MathUtils::seedRNG( const int&  a_seed )
{
   if(a_seed < 0) {
      std::random_device rd;
      global_rand_gen.seed(rd());
   }
   else {
      global_rand_gen.seed(a_seed);
   }
}

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

double MathUtils::gammainc( const Real  a_x,
                            const Real  a_a )
{
   // incomplete gamma function:  int_0^x(t^(a-1)*exp(-t^2))dt
   //
   
   // right now this is a Taylor expansion method assuming a=3/2
   // Make more general later
   //
   CH_assert(a_a==3./2.);

   Real soln = tgamma(a_a);
   if(a_x<10.0) {

      const int pmax = 41;
      Real coef;
      Real sign = -1.0;
      Real factorial = 1;

      soln = 0.0;
      for (auto p=1; p<pmax; p++) {
         coef = 0.5 + p;
         if(p>1) factorial = factorial*(p-1);
         sign = -sign;
       
         soln = soln + sign*pow(a_x,coef)/coef/factorial;
      }

   }

   return soln;

}

double MathUtils::rand()
{
   CH_TIME("MathUtils::rand()");
   
   static std::uniform_real_distribution<> dis(0,1);
   double randNum = dis(global_rand_gen);

   return randNum;

}

int MathUtils::randInt(const int& A, const int& B)
{
   CH_TIME("MathUtils::randInt()");
   
   //static std::random_device randDev;
   //static std::mt19937 twister(randDev());
   static std::uniform_int_distribution<int> dist;

   dist.param(std::uniform_int_distribution<int>::param_type(A,B));

   //return dist(twister);
   return dist(global_rand_gen);

}

double MathUtils::randn()
{
   CH_TIME("MathUtils::randn()");
   
   static std::normal_distribution<Real> disNorm(0.0,1.0);
   double randNum = disNorm(global_rand_gen);

   return randNum;

}

Real MathUtils::linearInterp( const std::vector<Real>&  a_X,
                              const std::vector<Real>&  a_Y,
                              const Real                a_X0,
                              const int                 a_index )
{
     
   Real Y0;
   CH_assert(a_index<a_X.size()-1);
       
   Y0 = ( a_Y[a_index+1]*(a_X0 - a_X[a_index]) 
      +   a_Y[a_index]*(a_X[a_index+1] - a_X0) )
      / ( a_X[a_index+1] - a_X[a_index] );

   return Y0;

}

Real MathUtils::linearInterp( int&                a_index,
                        const std::vector<Real>&  a_X,
                        const std::vector<Real>&  a_Y,
                        const Real                a_X0 )
{
     
   Real Y0;
   const int N = a_X.size();

   CH_assert(N==a_Y.size());        
   CH_assert(a_X0>=a_X.front());        
   CH_assert(a_X0<=a_X.back());        
 
   int index = N/2;
   while (a_X0 < a_X[index]) index--;
   while (a_X0 > a_X[index+1]) index++;
       
   Y0 = ( a_Y[index+1]*(a_X0 - a_Y[index]) 
      +   a_Y[index]*(a_X[index+1] - a_X0) )
      / ( a_X[index+1] - a_X[index] );

   a_index = index;
   return Y0;

}

#include "NamespaceFooter.H"
