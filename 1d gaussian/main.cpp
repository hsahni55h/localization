#include <iostream>
#include <math.h>

using namespace std;

double f(double mu, double sigma2, double x)
{
    //Use mu, sigma2 (sigma squared), and x to code the 1-dimensional Gaussian
    //Put your code here
    double prob = 1 / sqrt(2 * sigma2 * M_PI ) * exp( - 0.5 * pow((x-mu) , 2) / sigma2);
    return prob;
}

int main()
{
    cout << f(10.0, 4.0, 8.0) << endl;
    return 0;
}