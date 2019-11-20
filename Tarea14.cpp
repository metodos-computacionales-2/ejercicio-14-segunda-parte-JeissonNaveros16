#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

// Variable constantes globales
const double K = 100;
const double M = 2;
const double LAMBDA1 = 1;
const double LAMBDA3 = 3;
const double LAMBDA5 = 5;
const double DeltaT = 0.01;
const double GAMMA0 = 0; // constante para considerar la friccion
const double GAMMA1 = 1; // constante para considerar la friccion

// Declaracion de funciones
void writefile(string filename, double *t0, double *t1, double *t2, int n_points);

double f0(double t, double y0, double y1); // derivada de y0
double f1(double t, double y0, double y1, double lambda, double gamma); // derivada de y1
void rk4(double t, double h, double & y0, double & y1, double lambda, double gamma); // metodo de runge kutta 4 orden

double Euler0(double t, double y0, double y1);
double Euler1(double t, double y0, double y1, double lambda, double gamma);

int main(void)
{
    double *tiempo = new double[1000];
    double t = 0.0;
    for(int j = 0; j < 1000; j++)
    {
        *(tiempo + j) = t;
        t += DeltaT;
    }
    
    // Implementando Euler
    double *x0Euler = new double[1000];
    double *x1Euler = new double[1000];
    double t_0 = 0.0;
    double y_0 = 1.0;
    double y_1 = 0.0;
    cout << "\nPunto 3 - Solucion por metodo de Euler:" << endl;
    for(int i = 0; i < 1000; i++)
    {
        cout << t_0 << "\t" << y_0 << "\t" << y_1 << endl;
        t_0 += DeltaT;
        double temp0 = Euler0(t_0, y_0, y_1);
        double temp1 = Euler1(t_0, y_0, y_1, LAMBDA1, GAMMA0);
        *(x0Euler + i) = temp0;
        *(x1Euler + i) = temp1;
        y_0 = temp0;
        y_1 = temp1;
    }
    
    writefile("EulerSolution.dat", tiempo, x0Euler, x1Euler, 1000);
    
    // Implementando Runge-kutta
    double x, v, time;
    x = 1;
    v = 0;
    time = 0;
    double *x0RungeK = new double[1000];
    double *x1RungeK = new double[1000];
    cout << "\nPunto 5 - Solucion por metodo de Runge-kutta:" << endl;
    for(int i = 0; i < 1000; i++)
    {
        cout << time << "\t" << x << "\t" << v << endl;
        rk4(time, DeltaT, x, v, LAMBDA1, GAMMA0);
        *(x0RungeK + i) = x;
        *(x1RungeK + i) = v;
        time += DeltaT;
    }
    
    writefile("RungeKuttaSolution.dat", tiempo, x0RungeK, x1RungeK, 1000);
    
    // Ahora considerando la friccion
    double xf, vf, timef;
    xf = 1;
    vf = 0;
    timef = 0;
    double *x0RungeKFricc = new double[1000];
    double *x1RungeKFricc = new double[1000];
    cout << "\nPunto 7 - Solucion ecuacion diferencial con efecto de fricción:" << endl;
    for(int i = 0; i < 1000; i++)
    {
        cout << timef << "\t" << xf << "\t" << vf << endl;
        rk4(timef, DeltaT, xf, vf, LAMBDA1, GAMMA1);
        *(x0RungeKFricc + i) = xf;
        *(x1RungeKFricc + i) = vf;
        timef += DeltaT;
    }
    
    writefile("EDOFrictionSolution.dat", tiempo, x0RungeKFricc, x1RungeKFricc, 1000);
    
    // Ahora considerando un lambda de 3
    double xl3, vl3, time3;
    xl3 = 1;
    vl3 = 0;
    time3 = 0;
    double *x0l3 = new double[1000];
    double *x1l3 = new double[1000];
    cout << "\nPunto 5 - Solucion de EDO de Lambda=3:" << endl;
    for(int i = 0; i < 1000; i++)
    {
        cout << time3 << "\t" << xl3 << "\t" << vl3 << endl;
        rk4(time3, DeltaT, xl3, vl3, LAMBDA3, GAMMA0);
        *(x0l3 + i) = xl3;
        *(x1l3 + i) = vl3;
        time3 += DeltaT;
    }
    
    writefile("EDOLambda3Solution.dat", tiempo, x0l3, x1l3, 1000);
    
    // Ahora considerando un lambda de 5
    double xl5, vl5, time5;
    xl5 = 1;
    vl5 = 0;
    time5 = 0;
    double *x0l5 = new double[1000];
    double *x1l5 = new double[1000];
    cout << "\nPunto 8 - Solucion de EDO de Lambda=5:" << endl;
    for(int i = 0; i < 1000; i++)
    {
        cout << time5 << "\t" << xl5 << "\t" << vl5 << endl;
        rk4(time5, DeltaT, xl5, vl5, LAMBDA5, GAMMA0);
        *(x0l5 + i) = xl5;
        *(x1l5 + i) = vl5;
        time5 += DeltaT;
    }
    
    writefile("EDOLambda5Solution.dat", tiempo, x0l5, x1l5, 1000);
    
    delete[] x0Euler;
    delete[] x1Euler;
    delete[] x0RungeK;
    delete[] x1RungeK;
    delete[] x0RungeKFricc;
    delete[] x1RungeKFricc;
    delete[] x0l3;
    delete[] x1l3;
    delete[] x0l5;
    delete[] x1l5;
    
    return 0;
}

void writefile(string filename, double *t0, double *t1, double *t2, int n_points)
{
    // Aqui se abre el archivo con el nombre dado por parametro
    ofstream outfile;
    outfile.open(filename);
    
    // Luego se va imprimiendo en este los valores
    for (int j=0; j < n_points; j++)
    {
        outfile << *(t0 + j) << "\t" << *(t1 + j) << "\t" << *(t2 + j) << endl;
    }
    
    outfile.close();
}

double f0(double t, double y0, double y1)
{
    return y1;
}

double f1(double t, double y0, double y1, double lambda, double gamma)
{
    return (-K/M)*pow(y0, lambda) - (gamma/M)*y1;
}

void rk4(double t, double h, double & y0, double & y1, double lambda, double gamma) // metodo de runge kutta 4 orden
{
    double k10, k11, k20, k21, k30, k31, k40, k41;
    k10 = h*f0(t, y0, y1);
    k11 = h*f1(t, y0, y1, lambda, gamma);
    k20 = h*f0(t+h/2, y0 + k10/2, y1 + k11/2);
    k21 = h*f1(t+h/2, y0 + k10/2, y1 + k11/2, lambda, gamma);
    k30 = h*f0(t+h/2, y0 + k20/2, y1 + k21/2);
    k31 = h*f1(t+h/2, y0 + k20/2, y1 + k21/2, lambda, gamma);
    k40 = h*f0(t + h, y0 + k30, y1 + k31);
    k41 = h*f1(t + h, y0 + k30, y1 + k31, lambda, gamma);
    
    y0 = y0 + (1.0/6.0)*(k10 + 2*k20 + 2*k30 + k40);
    y1 = y1 + (1.0/6.0)*(k11 + 2*k21 + 2*k31 + k41);
}

// Metodos para resolver por Euler
double Euler0(double t, double y0, double y1)
{
    return y0 + DeltaT*f0(t, y0, y1);
}

double Euler1(double t, double y0, double y1, double lambda, double gamma)
{
    return y1 + DeltaT*f1(t, y0, y1, lambda, gamma);
}