#include "RecoLocalTracker/SiPixelRecHits/interface/VVIObj.cc"
#include<iostream>

double afun(x) { return x*x;}

int main() {
  double  sint=0;
  double  cint=0;
  for (double x=-20.5;x<20.5; x++2) {
    sincosint(x, sint, cint);
    std::cout << sint << " " << cint << " " 
	      << sinint(x) << " " << cosint(x) << std::endl;

    double x0, double rv;
    int res = dzero(-5,5,x0,rv,1.e-5,1000,afun);
    std::cout << res " " << x0 << " " << rv << " " std::endl;
    int res = dzero(5,-5,x0,rv,1.e-5,1000,afun);
    std::cout << res " " << x0 << " " << rv << " " std::endl;
    int res = dzero(5,15,x0,rv,1.e-5,1000,afun);
    std::cout << res " " << x0 << " " << rv << " " std::endl;

   return 0;

}
