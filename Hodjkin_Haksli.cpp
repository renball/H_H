#define _USE_MATH_DEFINES 
#include <cmath>
#include <math.h>
#include <iostream>

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

long double A_m(long double V) {
    //long double e = M_E;
    long double e;
    //double temp = (-1 * (V + 35) / 9);
    //e = pow(e, temp);
    e = exp(-1 * (V + 35 / 9));
    return (0.182 * (V + 35)) / (1 - e);
}

long double B_m(long double V) {
    long double e;// = M_E;
    //double temp = ((V + 35) / 9);
    //e = pow(e, temp);
    e = exp((V + 35) / 9);
    return (-0.124 * (V + 35)) / (1 - e);
}

long double A_n(long double V) {
    long double e; // = M_E;
    //double temp = (-1 * (V - 25) / 9);
    //e = pow(e,temp);
    e = exp(-1 * (V - 25) / 9);
    return (0.02 * (V - 25)) / (1 - e);
}

long double B_n(long double V) {
    long double e;// = M_E;
    //double temp = ((V - 25) / 9);
    //e = pow(e, temp);
    e = exp((V - 25) / 9);
    return (-0.002 * (V - 25)) / (1 - e);
}

long double A_h(long double V) {
    long double e;// = M_E;
    //double temp = (-1 * (V + 90) / 12);
    //e = pow(e, temp);
    e = exp(-1 * (V + 90) / 12);
    return 0.25 * e;
}

long double B_h(long double V) {
    //long double e = M_E;
    long double e1;
    long double e2;
    //double temp1 = ((V + 62) / 6);
    //double temp2 = ((V + 90) / 12);
    //e1 = pow(e, temp1);
    //e2 = pow(e, temp2);
    e1 = exp((V + 62) / 6);
    e2 = exp((V + 90) / 12);
    return 0.25 * (e1/e2);
}


long double dVdt(long double V, long  double m, long  double h, long  double n, long double I) {
    long double dV;
    dV = nu * (1 / C) * ((gNa * (m * m * m) * h * (ENa - V)) + (n*(gK * (Ek - V))) + (gLeak * (ELeak - V))+I);
    return dV;
}

long double dmdt(long double m, long double V) {
    long double am = A_m(V);
    long double bm = B_m(V);
    long double dm;
    dm = nu * (am * (1.0 - m) - bm * m);
    return dm;
}

long double dndt(long double n, long double V) {
    long double an = A_n(V);
    long double bn = B_n(V);
    long double dn;
    dn = nu * (an * (1.0 - n) - bn * n);
    return dn;
}

long double dhdt(long double h, long  double V) {
    long double ah = A_h(V);
    long double bh = B_h(V);
    long double dh;
    dh = nu * (ah * (1.0 - h) - bh * h);
    return dh;
}

//Runge_Kutta
long double RK_dvdt(long double V, long  double m, long  double h, long  double n, long double I) {
    long double yn = V;
    long double yn_1;
    long double k1;
    long double k2;
    long double k3;
    long double k4;
    long double hop=0.001;

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
    long double V[10000];
    V[0] = 14.8409;
    long double m = 0.9174;
    long double n = 0.0140;
    long double h = 0.0539;
    long double I = 0.9;
    long double t = 0.0;
    long double dt = 0.001;
    int i = 0;
    do  {
        cout << t << " "  << V[i] << " " << m << " " << n << " " << h << "\n";
        t += dt;
        /*
        V[i+1] = V[i] + dVdt(V[i], m, h, n, I) * dt;

        m = m+dmdt(m, V[i]) * dt;
        n = n+dndt(n, V[i]) * dt;
        h = h+dhdt(h, V[i]) * dt;
        */

        V[i + 1] = RK_dvdt(V[i], m, h, n, I);
        m = RK_dmdt(V[i], m);
        n = RK_dndt(V[i], n);
        h = RK_dhdt(V[i], h);

        i++;
    } while (t <= 1.0);

}
