#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

const double k = 0.7;
const double v0y =10;
const double v0x = 10;
const double angulo = 34.0;
const double g = 9.8;
const double dt = 0.01;
const double T = 2*v0y/g;

double Vy(double t);
double Vx();
double velocidady(double t, double y, double vy); 
double aceleraciony(double t, double y, double vy);
double velocidadx(double t, double x, double vx ); 
double aceleracionx(double t, double x, double vx); 
void rk4(double t, double dt, double &y, double &vy, double &x, double &vx);

int main()
{
	double x,vx,y,vy,time;
	time = 0.0;
	x = 0.0;
	y = 0.0;
	vx = Vx();
	vy = Vy(time);
	for(time = 0; time <= T; time += dt) 
  {
    cout <<" " << x << " " << y << endl;
    rk4(time, dt, y, vy, x, vx);
  }
  return 0;
}

double Vy(double t)
{
    return v0y - g*t;
}

double Vx()
{
    return v0x;
}

double velocidady(double t, double y0, double y1)
{
  return y1;
}

double aceleraciony(double t, double y0, double y1)
{
  return -g-k*Vy(t)*(Vy(t)/sqrt(Vx()*Vx() + Vy(t)*Vy(t)));
}

double velocidadx(double t,double x0,double x1)
{
    return x1;
}
double aceleracionx(double t,double x0, double x1)
{
    return -k*Vx()*(Vx()/sqrt(Vx()*Vx() + Vy(t)*Vy(t)));
}
void rk4(double t, double dt, double &y, double &vy, double &x, double &vx) 
{
    double k1y, k1vy, k2y, k2vy, k3y, k3vy, k4y, k4vy;
    double k1x, k1vx, k2x, k2vx, k3x, k3vx, k4x, k4vx;
    k1x = dt*velocidadx(t, x, vx);
    k1vx = dt*aceleracionx(t, x, vx);
    k2x = dt*velocidadx(t+dt/2, x + k1x/2, vx + k1vx/2);
    k2vx = dt*aceleracionx(t+dt/2, x + k1x/2, vx + k1vx/2);
    k3x = dt*velocidadx(t+dt/2, x + k2x/2, vx + k2vx/2);
    k3vx = dt*aceleracionx(t+dt/2, x + k2x/2, vx + k2vx/2);
    k4x = dt*velocidadx(t + dt,x + k3x, vx + k3vx);
    k4vx = dt*aceleracionx(t + dt, x + k3x, vx + k3vx);
    
    k1y = dt*velocidady(t, y, vy);
    k1vy = dt*aceleraciony(t, y, vy);
    k2y = dt*velocidady(t+dt/2, y + k1y/2, vy + k1vy/2);
    k2vy = dt*aceleraciony(t+dt/2, y + k1y/2, vy + k1vy/2);
    k3y = dt*velocidady(t+dt/2, y + k2y/2, vy + k2vy/2);
    k3vy = dt*aceleraciony(t+dt/2, y + k2y/2, vy + k2vy/2);
    k4y = dt*velocidady(t + dt, y + k3y, vy + k3vy);
    k4vy = dt*aceleraciony(t + dt, y + k3y, vy + k3vy);
    
    x = x + (1.0/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
	y = y + (1.0/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
    vx = vx + (1.0/6.0)*(k1vx + 2*k2vx + 2*k3vx + k4vx);
    vy = vy + (1.0/6.0)*(k1vy + 2*k2vy + 2*k3vy + k4vy);
}
