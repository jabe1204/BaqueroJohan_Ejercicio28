#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

const double k = 0.9;
const double masa = 10;
const double angulo = 35.0;
const double v0 = 22.0;
const double g = 9.8;
const double dt = 0.01;

double velocidady(double t, double y, double vy); 
double aceleraciony(double t, double y, double vy);
double velocidadx(double t, double x, double vx ); 
double aceleracionx(double t, double x, double vx); 
void rk4(double t, double dt, double &y, double &vy, double &x, double &vx);

int main()
{
  double x, v, time;
  x = 1;
  v = 0;
  for(time = 0; time <= 10; time += DeltaT) 
  {
    cout << time << "\t" << x << "\t" << v << endl;
    rk4(time, DeltaT, x, v);
  }
  return 0;
}

double v0y()
{
    return v0*sin(angulo*M_PI/180.0);
}

double Vely(double t)
{
    return v0y-g*t;
}

double Velx()
{
    return v0*cos(angulo*M_PI/180.0);
}

double velocidady(double t, double y0, double y1)
{
  return y1;
}

double aceleraciony(double t, double y0, double y1)
{
  return -g-k*Vely(t)*Vely(t)/sqrt(Velx*Velx + Vely(t)*Vely(t));
}

double velocidadx(double t,double x0,double x1)
{
    return x1;
}
double aceleracionx(double t,double x0, double x1)
{
    return -k*Velx*Velx/sqrt(Velx*Velx + Vely(t)*Vely(t));
}
void rk4(double t, double dt, double &y, double &vy, double &x, double &vx) 
{
  double k10, k11, k20, k21, k30, k31, k40, k41;
  k10 = dt*velocidady(t, y0, y1);
  k11 = dt*aceleraciony(t, y0, y1);
  k20 = dt*velocidady(t+h/2, y0 + k10/2, y1 + k11/2);
  k21 = dt*aceleraciony(t+h/2, y0 + k10/2, y1 + k11/2);
  k30 = dt*velocidady(t+h/2, y0 + k20/2, y1 + k21/2);
  k31 = dt*aceleraciony(t+h/2, y0 + k20/2, y1 + k21/2);
  k40 = dt*velocidady(t + h, y0 + k30, y1 + k31);
  k41 = dt*aceleraciony(t + h, y0 + k30, y1 + k31);

  y = y + (1.0/6.0)*(k10 + 2*k20 + 2*k30 + k40);
  vy = vy + (1.0/6.0)*(k11 + 2*k21 + 2*k31 + k41);
}