#include <iostream>
#include <cstdlib>
#include <random>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

mt19937 gen(547426715462);
uniform_real_distribution<double> ran_u(0.0, 1.0);
normal_distribution<double> ran_g(0.0, 1.0);

//Constructs the network given the number of patches
void construct_network(vector< vector<double> > &network, const double M)
{
    int i,j;

    network = vector< vector<double> >(M, vector<double>(M));

    double x1, y1, x2, y2;
    vector<double> x(M);
    vector<double> y(M);

    for (i=0; i < M; i++)
    {
        x[i] = ran_u(gen);
        y[i] = ran_u(gen);
    }

    for (i=0; i < M; i++)
    {
        network[i][i] = 0.0;
        for (j=i+1; j < M; j++)
        {
            network[i][j] = 1.0/((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]));
            network[j][i] = 1.0/((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]));
        }
    }

    return;
}

double kernel(const double x, const double thres)
{
    return (x < thres) / thres;
}

double sigmoidal(const double x)
{
    return 1.0 / (1.0 + exp(x));
}

void initial_conditions(const int M, vector<double> &rho,  vector<double> &R, const double k0, vector<double> &k, vector<double> &h)
{
    int i;

    rho = vector<double>(M, 0.0);
    R = vector<double>(M);

    k = vector<double>(M);
    h = vector<double>(M, 1.0);

    for (i=0; i < M; i++)
    {
        R[i] = 1.0;
        k[i] = k0 + abs(ran_g(gen)); 
    }

    rho[0] = 0.5*R[0];

    return;    
}

void simulation_step(const int M, const vector< vector<double> > &network, vector<double> &rho, vector<double> &R, 
                    const double d, const double lambda, const double sigma, const double theta, 
                    const vector<double> &k, const vector<double> h, const double dt, const double sqdt)
{
    int i,j;
    double input, aux_input;
    double resource_sum;

    double deterministic;
    double noise;

    vector<double> old_rho(M);
    vector<double> old_R(M);
    swap(old_R, R);
    swap(old_rho, rho);

    //Main loop over patches 
    for (i=0; i < M; i++)
    {
        //Input from neighbouring patches 
        input = 0.0;
        for (j=0; j < M; j++)
        {
            aux_input = lambda * sigmoidal(k[j]);
            input += network[i][j] * aux_input * old_rho[j] * kernel(old_rho[j], theta);  
            //resource_sum += k[j]*old_rho[j]*old_R[j];
        }
        input *= 1.0 / network[i].size();
        //input = 0.0;

        //Deterministic part for the density
        deterministic = (k[i]*old_rho[i] + input) * (1.0 - old_rho[i]/old_R[i]) - d*old_rho[i];

        rho[i] = old_rho[i] + dt*deterministic + sqdt*sqrt(old_rho[i])*sigma*ran_g(gen); 
        rho[i] = max(rho[i], 0.0);                   
        R[i] = old_R[i] + dt * (h[i] * old_R[i] * (1.0 - old_R[i])  - k[i]*old_rho[i]*old_R[i]);
        //R[i] = old_R[i] + dt * (h[i] - k[i]*old_rho[i]*old_R[i]);
    }

    return;
}

int main(int argc, char* argv[])
{
    int i;

    int M = 100; //Number of patches

    double t;
    double tf = 100.0;
    const double dt = 0.01;
    const double sqdt = sqrt(dt);

    ofstream output;
    int writetime, ctime;
    writetime = 10;
    string file;


    vector< vector<double> > network;
    vector<double> rho, R;
    vector<double> h, k; 
    double k0, lambda, theta, d, sigma;

    if (argc == 7)
    {
        k0     = stod(argv[1]);
        lambda = stod(argv[2]);
        d      = stod(argv[3]);
        theta  = stod(argv[4]);
        sigma  = stod(argv[5]);
        file   = string(argv[6]);
    }
    else 
    {
        cout << "Incorrect number of parameters" << endl;
        return EXIT_SUCCESS;
    }

    construct_network(network, M);
    initial_conditions(M, rho, R, k0, k, h);

    output.open(file);
    ctime = 0;
    for (t=0.0; t < tf; t += dt)
    {
        simulation_step(M, network, rho, R, d, lambda, sigma, theta, k, h, dt, sqdt);
        if (ctime % writetime == 0)
        {
            for (i=0; i < M; i++) output << rho[i] << " " << R[i] << " ";
            output << endl;
        }
        ctime++;
    }
    output.close();

    return EXIT_SUCCESS;
}