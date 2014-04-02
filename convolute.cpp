#include <iostream>
#include <vector>
 
#include <cmath>

using namespace std;

double Gauss(const double &t, const double &s) {
    double amp = 1.0;
    double mu = 0.0;
    double offset = 0.0;

    double coeff = amp / (s*sqrt(2*M_PI));
    double exponent = exp(-pow(t-mu,2) / 2 / s / s); 
    return(coeff * exponent + offset);
}

double GaussStarGauss(const double &t, const double &tau, const double &s1, 
                      const double &s2) {
    return(Gauss(tau,s1)*Gauss(t-tau,s2));
}

double AdaptiveSimpsonsAux(const double &a, const double &b, 
                           const double &epsilon, const double &S, 
                           const double &fa, const double &fb, 
                           const double &fc, const int &bottom, 
                           const double &t, const double &s1, 
                           const double &s2) {
    double c = (a + b)/2, h = b - a;
    double d = (a + c)/2, e = (c + b)/2;
    double fd = GaussStarGauss(t,d,s1,s2), fe = GaussStarGauss(t,e,s1,s2);
    double Sleft = (h/12)*(fa + 4*fd + fc);
    double Sright = (h/12)*(fc + 4*fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)
        return( S2 + (S2 - S)/15 );
    return( AdaptiveSimpsonsAux(a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1, 
                                t,s1,s2) + 
            AdaptiveSimpsonsAux(c, b, epsilon/2, Sright, fc, fb, fe, bottom-1, 
                                t,s1,s2) );
}

double AdaptiveSimpsons(const double &a, const double &b, // interval 
                        const double &epsilon, // error tolerance
                        const int &maxRecursionDepth, // recursion cap
                        const double &t, const double &s1, const double &s2){ 
    double c = (a + b)/2, h = b - a;
    double fa = GaussStarGauss(t,a,s1,s2), fb = GaussStarGauss(t,b,s1,s2), 
        fc = GaussStarGauss(t,c,s1,s2);
    double S = (h/6)*(fa + 4*fc + fb);
    return AdaptiveSimpsonsAux(a, b, epsilon, S, fa, fb, fc, maxRecursionDepth, 
                               t, s1, s2);
}

int main(int argc, char* argv[]) {
    for(double t = -10; t < 10; t+=0.1) {
        cout << t << " " << Gauss(t, 1.0) << " " << Gauss(t, 2.0) << " ";
        cout << AdaptiveSimpsons(-1e6, 1e6, 1e-18, 30, t, 1.0, 2.0) << endl;
    } //for (double t
}
