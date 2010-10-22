#include "RecoLocalTracker/SiPixelRecHits/src/VVIObj.cc"
#include<iostream>

double afun(double x) { return x;}

int main() {
  using namespace VVIObjDetails;

  double  sint=0;
  double  cint=0;
  for (double x=-20.5;x<20.5; x+=2) {
    sincosint(x, sint, cint);
    std::cout << sint << " " << cint << " " 
	      << sinint(x) << " " << cosint(x) << std::endl;
  }
    double x0; double rv;
    int res = dzero(-5,5,x0,rv,1.e-5,1000,afun);
    std::cout << res << " " << x0 << " " << rv << " " << std::endl;
    res = dzero(5,-5,x0,rv,1.e-5,1000,afun);
    std::cout << res << " " << x0 << " " << rv << " " << std::endl;
    res = dzero(5,15,x0,rv,1.e-5,1000,afun);
    std::cout << res << " " << x0 << " " << rv << " " << std::endl;

    double kappa[] = { 0.0681354, 0.0725822, 2.27383  };
    double xvav[]  = { 0.434149,  7.45804, 2.27383  };

    double beta2=1;
    for (int i=0; i!=3; ++i) {
      VVIObj vvidist(kappa[i], beta2, 1);
      double prvav = vvidist.fcn(xvav[i]);
      std::cout << "vav: " << kappa << " " << xvav << " " << prvav
               << " " << TMath::VavilovI(xvav[i], kappa[i], beta2) << std::endl;
    }

   return 0;

}
