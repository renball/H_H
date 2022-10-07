#define _USE_MATH_DEFINES 
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

const int nu = 1000;
const int C = 1;
const int gNa = 40;
const int gK = 35;
const float gLeak = 0.3;
const int ENa = 55;
const int Ek = -77;
const float ELeak = -54.4;
const float ksyn = 0.2;

double A_m(double V) {
    return (0.182 * (V + 35)) / (1 - exp(-1 * (V + 35 / 9)));
}

double B_m(double V) {
    return (-0.124 * (V + 35)) / (1 - exp((V + 35) / 9));
}

double A_n(double V) {
    return (0.02 * (V - 25)) / (1 - exp(-1 * (V - 25) / 9));
}

double B_n(double V) {
    return (-0.002 * (V - 25)) / (1 - exp((V - 25) / 9));
}

double A_h(double V) {
    return 0.25 * exp(-1 * (V + 90) / 12);
}

double B_h(double V) {
    return 0.25 * (exp((V + 62) / 6) / exp((V + 90) / 12));
}


double dVdt(double V, double m, double h, double n, double I) {

    return nu * (1 / C) * ((gNa * (m * m * m) * h * (ENa - V)) + (n * (gK * (Ek - V))) + (gLeak * (ELeak - V)) + I);
}

double dmdt(double m, double V) {
    double am = A_m(V);
    double bm = B_m(V);
    return nu * (am * (1.0 - m) - bm * m);;
}

double dndt(double n, double V) {
    double an = A_n(V);
    double bn = B_n(V);
    return nu * (an * (1.0 - n) - bn * n);
}

double dhdt(double h, double V) {
    double ah = A_h(V);
    double bh = B_h(V);
    return nu * (ah * (1.0 - h) - bh * h);
}

//Runge_Kutta
long double RK_dvdt(long double V, long  double m, long  double h, long  double n, long double I) {
    long double yn = V;
    long double yn_1;
    long double k1;
    long double k2;
    long double k3;
    long double k4;
    long double hop = 0.001;

    k1 = dVdt(yn, m, h, n, I);
    k2 = dVdt((yn + hop * (k1 / 2)), m, h, n, I);
    k3 = dVdt((yn + hop * (k2 / 2)), m, h, n, I);
    k4 = dVdt((yn + hop * k3), m, h, n, I);

    return yn + (k1 + 2 * k2 + 2 * k3 + k4) * hop / 6;
}

long double RK_dndt(long double V, long  double n) {
    long double yn = n;
    long double yn_1;
    long double k1;
    long double k2;
    long double k3;
    long double k4;
    long double hop = 0.001;

    k1 = dndt(yn, V);
    k2 = dndt((yn + hop * (k1 / 2)), V);
    k3 = dndt((yn + hop * (k2 / 2)), V);
    k4 = dndt((yn + hop * k3), V);

    return yn + (k1 + 2 * k2 + 2 * k3 + k4) * hop / 6;
}

long double RK_dmdt(long double V, long  double m) {
    long double yn = m;
    long double yn_1;
    long double k1;
    long double k2;
    long double k3;
    long double k4;
    long double hop = 0.001;

    k1 = dndt(yn, V);
    k2 = dndt((yn + hop * (k1 / 2)), V);
    k3 = dndt((yn + hop * (k2 / 2)), V);
    k4 = dndt((yn + hop * k3), V);

    return yn + (k1 + 2 * k2 + 2 * k3 + k4) * hop / 6;
}

long double RK_dhdt(long double V, long  double h) {
    long double yn = h;
    long double yn_1;
    long double k1;
    long double k2;
    long double k3;
    long double k4;
    long double hop = 0.001;

    k1 = dndt(yn, V);
    k2 = dndt((yn + hop * (k1 / 2)), V);
    k3 = dndt((yn + hop * (k2 / 2)), V);
    k4 = dndt((yn + hop * k3), V);

    return yn + (k1 + 2 * k2 + 2 * k3 + k4) * hop / 6;
}

int main()
{
    fstream tx;
    tx.open("Results.txt");

    double VC[4]; //Current
    double VN[4]; //Next
    VC[0] = 14.8409; //V
    VC[1] = 0.9174; //m
    VC[2] = 0.0140; //n
    VC[3] = 0.0539; //h
    double I = 0.9;
    double t = 0.0;
    double dt = 0.00005;

    while (t <= 1.0) {
        //cout << t << " "  << VC[0] << " " << VC[1] << " " << VC[2] << " " << VC[3] << "\n";
        tx << t << "     " << VC[0] << "     " << VC[1] << "     " << VC[2] << "     " << VC[3] << endl;

        VN[0] = VC[0] + dVdt(VC[0], VC[1], VC[3], VC[2], I) * dt;

        VN[1] = VC[1] + dmdt(VC[1], VC[0]) * dt;
        VN[2] = VC[2] + dndt(VC[2], VC[0]) * dt;
        VN[3] = VC[3] + dhdt(VC[3], VC[0]) * dt;

        VC[0] = VN[0];
        VC[1] = VN[1];
        VC[2] = VN[2];
        VC[3] = VN[3];

        t += dt;
        /*
        V[i + 1] = RK_dvdt(V[i], m, h, n, I);
        m = RK_dmdt(V[i], m);
        n = RK_dndt(V[i], n);
        h = RK_dhdt(V[i], h);
        */
    }
    tx.close();

}
