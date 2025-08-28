#ifndef MAIN_H
#define MAIN_H

typedef struct {
    double x, y, z, r;
    double vx, vy, vz;
    double t;
    int bounces;
    double magtron_phase;
} state;

double AzimuthalAngle(double, double);
void BounceChecker(double, state*);
double CyclotronRadius(double, double, double);
double dvx_dt (double, double, double, double);
double dvy_dt (double, double, double, double);
double dvz_dt (double, double, double, double);
double dx_dt (double);
double dy_dt (double);
double dz_dt (double);
void MagnetronMotion(double, state*);
void PositionRandomizer(double*, double*, double*, double*);
double RangedRandomizer(double, double);
state RungeKutta (state, double, double, double);
void SaveParameters(state*, state*);
double ScalarPotentialApproximation(double, double, double);
double TotalEnergy (double, double, double, double, double, double);
void VectorPotentialCalculation(double, double, double, double*, double*, double*);
void VelocityRandomizer(double*, double*, double*, double);

#endif