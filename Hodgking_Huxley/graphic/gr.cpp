#define _USE_MATH_DEFINES 
#include <cmath>
#include <math.h>
#include "gr.h"
#include <Windows.h>
using namespace graphic;

int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	Application::Run(gcnew gr);
	return 0;
}

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


double dVdt(double V, double m, double n, double h, double I) {

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


System::Void graphic::gr::button1_Click(System::Object^ sender, System::EventArgs^ e)
{

	I = Convert::ToDouble(textBox1->Text);

	this->chart1->Series[0]->Points->Clear();

    double VC[4]; //Current
    double VN[4]; //Next
    VC[0] = 14.8409; //V
    VC[1] = 0.9174; //m
    VC[2] = 0.0140; //n
    VC[3] = 0.0539; //h
    //double I = 0.9;
    double t = 0.0;
    double dt = 0.00005;

    while (t <= 0.05) {
       
        this->chart1->Series[0]->Points->AddXY(t, VC[0]);

        VN[0] = VC[0] + dVdt(VC[0], VC[1], VC[2], VC[3], I) * dt;

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

	return System::Void();
}

System::Void graphic::gr::chart1_Click(System::Object^ sender, System::EventArgs^ e)
{
	return System::Void();
}

