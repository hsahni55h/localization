#include <iostream>

using namespace std;

double w[] = { 0.6, 1.2, 2.4, 0.6, 1.2 };

//Define a  ComputeProb function and compute the Probabilities

double sum = 0;
int n = sizeof(w)/sizeof(w[0]);

void ComputeProb(double w[], int n)
{
    for(int i = 0; i < n; i++)
    {
        sum = sum + w[i];
    }

    for(int j = 0; j < n; j++)
    {
        w[j] =  w[j]/sum;
        cout << "probabilities of particle" << j + 1 << "=" << w[j] << endl;
    }
}


int main()
{
    ComputeProb(w, n);
    return 0;
}


