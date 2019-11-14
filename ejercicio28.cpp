#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
/*Basado en https://josemontanac.github.io/Laboratorio-Metodos-Computacionales/2/Ordinary%20differential%20equations%20(ODE).slides.html#/10*/

const double k = 0.7;
const double v0y =10;
const double v0x = 10;
const double angulo = 34.0;
const double g = 9.8;
const double dt = 0.01;
const double T = 2*v0y/g;

double velocidady(double t,double y, double vy, double x,double vx); 
double aceleraciony(double t,double y, double vy, double x,double vx);
double velocidadx(double t,double y, double vy, double x,double vx); 
double aceleracionx(double t,double y, double vy, double x,double vx); 
void rk4(double t, double dt, double &y, double &vy, double &x, double &vx);

int main()
{
	double x,vx,y,vy,time,xs,ys;
	time = 0.0;
	x = 0.0;
	y = 0.0;
	vx = v0y;
	vy = v0y;
	
	for(time = 0; time <= T; time += dt) 
  {
		xs = v0x*time;
		ys = v0y*time - g*time*time/2;
		cout <<" " << x << " " << y << " "<<xs<<" "<<ys<< endl;
		rk4(time, dt, y, vy, x, vx);
  }
  return 0;
}

double velocidady(double t, double y, double vy, double x, double vx)
{
  return vy;
}

double aceleraciony(double t, double y, double vy,double x, double vx)
{
  return -g-k*vy*vy/sqrt(vx*vx + vy*vy);
}

double velocidadx(double t,double y, double vy, double x,double vx)
{
    return vx;
}
double aceleracionx(double t,double y, double vy, double x, double vx)
{
    return -k*vx*vx/sqrt(vx*vx + vy*vy);
}
void rk4(double t, double dt, double &y, double &vy, double &x, double &vx) 
{
    double k1y, k1vy, k2y, k2vy, k3y, k3vy, k4y, k4vy;
    double k1x, k1vx, k2x, k2vx, k3x, k3vx, k4x, k4vx;
    k1x = dt*velocidadx(t, y, vy, x, vx);
    k1vx = dt*aceleracionx(t,y, vy, x, vx);
    k2x = dt*velocidadx(t+dt/2,y, vy, x + k1x/2, vx + k1vx/2);
    k2vx = dt*aceleracionx(t+dt/2,y, vy, x + k1x/2, vx + k1vx/2);
    k3x = dt*velocidadx(t+dt/2,y, vy, x + k2x/2, vx + k2vx/2);
    k3vx = dt*aceleracionx(t+dt/2,y, vy, x + k2x/2, vx + k2vx/2);
    k4x = dt*velocidadx(t + dt,y, vy,x + k3x, vx + k3vx);
    k4vx = dt*aceleracionx(t + dt,y, vy, x + k3x, vx + k3vx);
    
    k1y = dt*velocidady(t, y, vy, x, vx);
    k1vy = dt*aceleraciony(t, y, vy, x, vx);
    k2y = dt*velocidady(t+dt/2, y + k1y/2, vy + k1vy/2, x, vx);
    k2vy = dt*aceleraciony(t+dt/2, y + k1y/2, vy + k1vy/2, x, vx);
    k3y = dt*velocidady(t+dt/2, y + k2y/2, vy + k2vy/2, x, vx);
    k3vy = dt*aceleraciony(t+dt/2, y + k2y/2, vy + k2vy/2, x, vx);
    k4y = dt*velocidady(t + dt, y + k3y, vy + k3vy, x, vx);
    k4vy = dt*aceleraciony(t + dt, y + k3y, vy + k3vy, x, vx);
    
    x = x + (1.0/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
	y = y + (1.0/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
    vx = vx + (1.0/6.0)*(k1vx + 2*k2vx + 2*k3vx + k4vx);
    vy = vy + (1.0/6.0)*(k1vy + 2*k2vy + 2*k3vy + k4vy);
}
