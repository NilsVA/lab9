#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N);
void initialize(double* const u, const double dx, const double xmin,
                const int N);
void step(double* u0, double* u1, const int N, const double V, const double dx, double dt);
//---------------------------------------
int main(){

  const double tEnd = 5.;
  const double V = 1.;

  const int N  = 256;
  const double xmin = -10;
  const double xmax =  10;
  const double dx = (xmax-xmin)/(N-1);
  double dt = 0.5*dx/V;
  const int Na = 10; // Number of output files up to tEnd
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  double* h;

  stringstream strm;

  initialize(u0,dx, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N);

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

      step(u0,u1,N,V,dx,dt); 	// Put call to step function here

      h=u1;
      u1=u0;// swap arrays u0 <-> u1,
      u0=h;  // however do not copy values, be more clever ;)
   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N);
  }


  delete[] u0;
  delete[] u1;
  return 0;
}
//-----------------------------------------------

//-----------------------------------------------
void initialize(double* const u, const double dx, const double xmin,
                const int N)
{
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     if (fabs(x)<=1.0)
       u[i] = 1;
     else
      u[i] =0;
   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] << endl;
   }
   out.close();
}

void step(double* u0, double* u1, const int N, const double V, const double dx, double dt)
{
   u1[0]= (-1.)*(dt*V/(2*dx)) * (u0[1]) + u0[0];
  for(int i=1; i<N-1; i++)
   {
     u1[i]= (-1.)*(dt*V/(2*dx)) * (u0[i+1] - u0[i-1]) + u0[i];
     
   }
  u1[N-1]= (-1.)*(dt*V/(2*dx)) * (u0[N-2]) + u0[N-1];
  
}