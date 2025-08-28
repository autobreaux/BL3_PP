#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "main.h"
#include "constants.h"

// To Do:
// Calculate z-cor of max potential using cubic spline
// Randomize KE values using Vonemuann Sampling
// Make array of scalar potential values
// Complete parameter subroutines and saving process
    // Think about making two outfiles: one for rejected protons and one for trapper protons

// List of parameters to save:
// Starting cyclotron radius, final cyclotron radius, starting azimuthal angle, final azimuthal angle, number of z bounces,
// number of magnetron orbits, total time, etc.

///////////////////
// Main Function //
///////////////////

int main()
{
    // Initial Value Declarations //
    state initial;
    double initial_KE = 450; // Units of eV (This will be randomized later) //
    
    // Randomized Position & Velocity Initialization //
    srand(time(NULL));
    PositionRandomizer(&initial.x, &initial.y, &initial.z, &initial.r);
    VelocityRandomizer(&initial.vx, &initial.vy, &initial.vz, initial_KE);

    // Fluctuating Value Declarations //
    state proton = initial;
    double Ex = Ex0, Ey = Ey0, Ez = Ez0;

    // Run Main Runge Kutta //
    proton = RungeKutta(proton, Ex, Ey, Ez);

    // Write Parameters into File //
    SaveParameters(&initial, &proton);

    return 0;
}

/////////////////////////////
// Runge Kutta Subroutines //
/////////////////////////////

    // Differential Equation Subroutines //
    double dvx_dt (double vx, double vy, double vz, double Ex)
    {
        return ((q/mp) * (Ex + (vy * Bz) - (vz * By)));
    }

    double dvy_dt (double vx, double vy, double vz, double Ey)
    {
        return ((q/mp) * (Ey + (vz * Bx) - (vx * Bz)));
    }

    double dvz_dt (double vx, double vy, double vz, double Ez)
    {
        return ((q/mp) * (Ez + (vx * By) - (vy * Bx)));
    }

    double dx_dt (double vx)
    {
        return vx;
    }

    double dy_dt (double vy)
    {
        return vy;
    }

    double dz_dt (double vz)
    {
        return vz;
    }

    // Main Runge Kutta Subroutine //
    state RungeKutta (state proton, double Ex, double Ey, double Ez)
    {
        // Calculate Scalar Potentials of Boundary Walls //
        double maxpotential = ScalarPotentialApproximation (proton.x, proton.y, z1_maxpotential) * q;

        // Main Loop //
        for (int i = 0; i < n; i++)
        {
            // Placeholder Variables (Bounces, Magnetron Motion) //
            double vz_placeholder = proton.vz;
            double theta_placeholder = AzimuthalAngle(proton.x, proton.y);

            // MAKE THIS A SUBROUTINE
            // Check Total Energy and Bounds //
            if (proton.z <= z1_maxpotential || proton.z >= z2_maxpotential)
            {
                double TotE = TotalEnergy(proton.x, proton.y, proton.z, proton.vx, proton.vy, proton.vz);
                if (TotE < maxpotential)
                {
                    proton.vz = -proton.vz;
                }
            }

            // Retrieve Electric Field Vector Components //
            VectorPotentialCalculation(proton.x, proton.y, proton.z, &Ex, &Ey, &Ez);

            // Perform Stepping //
            double k1_vx = h * dvx_dt(proton.vx, proton.vy, proton.vz, Ex);
            double k1_vy = h * dvy_dt(proton.vx, proton.vy, proton.vz, Ey);
            double k1_vz = h * dvz_dt(proton.vx, proton.vy, proton.vz, Ez);
            double k1_x = h * dx_dt(proton.vx);
            double k1_y = h * dy_dt(proton.vy);
            double k1_z = h * dz_dt(proton.vz);

            double k2_vx = h * dvx_dt(proton.vx + k1_vx / 2, proton.vy + k1_vy / 2, proton.vz + k1_vz / 2, Ex);
            double k2_vy = h * dvy_dt(proton.vx + k1_vx / 2, proton.vy + k1_vy / 2, proton.vz + k1_vz / 2, Ey);
            double k2_vz = h * dvz_dt(proton.vx + k1_vx / 2, proton.vy + k1_vy / 2, proton.vz + k1_vz / 2, Ez);
            double k2_x = h * dx_dt(proton.vx + k1_vx / 2);
            double k2_y = h * dy_dt(proton.vy + k1_vy / 2);
            double k2_z = h * dz_dt(proton.vz + k1_vz / 2);

            double k3_vx = h * dvx_dt(proton.vx + k2_vx / 2, proton.vy + k2_vy / 2, proton.vz + k2_vz / 2, Ex);
            double k3_vy = h * dvy_dt(proton.vx + k2_vx / 2, proton.vy + k2_vy / 2, proton.vz + k2_vz / 2, Ey);
            double k3_vz = h * dvz_dt(proton.vx + k2_vx / 2, proton.vy + k2_vy / 2, proton.vz + k2_vz / 2, Ez);
            double k3_x = h * dx_dt(proton.vx + k2_vx / 2);
            double k3_y = h * dy_dt(proton.vy + k2_vy / 2);
            double k3_z = h * dz_dt(proton.vz + k2_vz / 2);

            double k4_vx = h * dvx_dt(proton.vx + k3_vx, proton.vy + k3_vy, proton.vz + k3_vz, Ex);
            double k4_vy = h * dvy_dt(proton.vx + k3_vx, proton.vy + k3_vy, proton.vz + k3_vz, Ey);
            double k4_vz = h * dvz_dt(proton.vx + k3_vx, proton.vy + k3_vy, proton.vz + k3_vz, Ez);
            double k4_x = h * dx_dt(proton.vx + k3_vx);
            double k4_y = h * dy_dt(proton.vy + k3_vy);
            double k4_z = h * dz_dt(proton.vz + k3_vz);

            proton.vx += (1.0/6.0) * (k1_vx + (2 * k2_vx) + (2 * k3_vx) + k4_vx);
            proton.vy += (1.0/6.0) * (k1_vy + (2 * k2_vy) + (2 * k3_vy) + k4_vy);
            proton.vz += (1.0/6.0) * (k1_vz + (2 * k2_vz) + (2 * k3_vz) + k4_vz);
            proton.x += (1.0/6.0) * (k1_x + (2 * k2_x) + (2 * k3_x) + k4_x);
            proton.y += (1.0/6.0) * (k1_y + (2 * k2_y) + (2 * k3_y) + k4_y);
            proton.z += (1.0/6.0) * (k1_z + (2 * k2_z) + (2 * k3_z) + k4_z);
            proton.t += h;
            proton.r = proton.x * proton.x + proton.y * proton.y;

            // Check for Bounce //
            BounceChecker(vz_placeholder, &proton);

            // Track Magnetron Motion //
            MagnetronMotion(theta_placeholder, &proton);
        }
        return proton;
    }

//////////////////////////////////////////////
// Electric Potential and Field Subroutines //
//////////////////////////////////////////////

    // Scalar Electric Potential Approximation Subroutine //
    double ScalarPotentialApproximation(double x, double y, double z)
    {
        // Cylinder Voltages in Order from Left to Right //
        double voltages[] = {900.0, 900.0, 900.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 900.0, 900.0, 900.0};
        int size_of_voltages = sizeof(voltages)/sizeof(voltages[0]);

        double midpoint_affwell = L / 2.0;
        double midpoint_cancwell = (L / 2.0) * -1.0;

        double potential_affwell = 0.0;
        double potential_cancwell = 0.0;
        double potential = 0.0;

        // Primary Loop for Subroutine //
        for(int i = 0; i < size_of_voltages; i++)
        {
            // Definitions for the Approximation //
            double scalar_v1_affwell = voltages[i];
            double scalar_v2_affwell = 0.0;

            double scalar_v1_cancwell = -voltages[i];
            double scalar_v2_cancwell = 0.0;

            // First Considering the Affecting Well //
            // if (fabs((omega / R) * (z - midpoint_affwell)) <= 1.0)
            // {
            potential_affwell = (0.5) * (scalar_v1_affwell + scalar_v2_affwell) + (0.5) * (scalar_v2_affwell - scalar_v1_affwell) * ((((2.0/3.0) * (1.0 / (2.0 * omega * (L / (2.0 * R)))) * log((cosh(omega * (((z - midpoint_affwell) + (L/2.0)) / R)))/(cosh(omega * (((z - midpoint_affwell) - (L/2.0)) / R))))) + ((1.0/3.0) * (4.0 / (M_PI * M_PI * omega * ((L / 2.0) / R))) * atan((sinh((M_PI * omega * (z - midpoint_affwell)) / (2.0 * R))) / (cosh((M_PI * omega * (L / 2.0)) / (2.0 * R)))) * atan((sinh((M_PI * omega * (L / 2.0)) / (2.0 * R))) / (cosh((M_PI * omega * (z - midpoint_affwell)) / (2.0 * R)))) * cosh((M_PI * omega * (z - midpoint_affwell)) / (2.0 * R)))) + ((((x * x + y * y) * omega) / (2.0 * R * R)) * (1.0 / cosh((omega * (z - midpoint_affwell)) / R)) * (1.0 / cosh((omega * (z - midpoint_affwell)) / R)) * tanh((omega * (z - midpoint_affwell)) / R)));
            // }
            // if (fabs((omega / R) * (z - midpoint_affwell)) > 1.0)
            // {
            //     if (z - midpoint_affwell > 0)
            //     {
            //         potential_affwell = (0.5) * (scalar_v1_affwell + scalar_v2_affwell) + (0.5) * (scalar_v2_affwell - scalar_v1_affwell) * ((((2.0/3.0) * (1.0 / (2.0 * omega * (L / (2.0 * R)))) * log((cosh(omega * (((z - midpoint_affwell) + (L/2.0)) / R)))/(cosh(omega * (((z - midpoint_affwell) - (L/2.0)) / R))))) + ((1.0/3.0) * (1.0 - ((4.0 / (M_PI * M_PI * omega * ((L / 2.0) / R))) * atan(1.0 / (sinh((M_PI * omega * (z - midpoint_affwell)) / (2.0 * R)))) * atan((sinh((M_PI * omega * (L / 2.0)) / (2.0 * R))) / (cosh((M_PI * omega * (z - midpoint_affwell)) / (2.0 * R)))) * cosh((M_PI * omega * (z - midpoint_affwell)) / (2.0 * R)))))) + ((((x * x + y * y) * omega) / (2.0 * R * R)) * (1.0 / cosh((omega * (z - midpoint_affwell)) / R)) * (1.0 / cosh((omega * (z - midpoint_affwell)) / R)) * tanh((omega * (z - midpoint_affwell)) / R)));
            //     }
            //     if (z - midpoint_affwell < 0)
            //     {
            //         potential_affwell = (0.5) * (scalar_v1_affwell + scalar_v2_affwell) + (0.5) * (scalar_v2_affwell + scalar_v1_affwell) * ((((2.0/3.0) * (1.0 / (2.0 * omega * (L / (2.0 * R)))) * log((cosh(omega * (((-1.0 * (z - midpoint_affwell)) + (L/2.0)) / R)))/(cosh(omega * (((-1.0 * (z - midpoint_affwell)) - (L/2.0)) / R))))) + ((1.0/3.0) * (1.0 - ((4.0 / (M_PI * M_PI * omega * ((L / 2.0) / R))) * atan(1.0 / (sinh((M_PI * omega * (-1.0 * (z - midpoint_affwell))) / (2.0 * R)))) * atan((sinh((M_PI * omega * (L / 2.0)) / (2.0 * R))) / (cosh((M_PI * omega * (-1.0 * (z - midpoint_affwell))) / (2.0 * R)))) * cosh((M_PI * omega * (-1.0 * (z - midpoint_affwell))) / (2.0 * R)))))) - ((((x * x + y * y) * omega) / (2.0 * R * R)) * (1.0 / cosh((omega * (-1.0 * (z - midpoint_affwell))) / R)) * (1.0 / cosh((omega * (-1.0 * (z - midpoint_affwell))) / R)) * tanh((omega * (-1.0 * (z - midpoint_affwell))) / R)));
            //     }
            // }

            // // Second Consider the Cancelling Well //
            // if (fabs((omega / R) * (z - midpoint_cancwell)) <= 1.0)
            // {
            potential_cancwell = (0.5) * (scalar_v1_cancwell + scalar_v2_cancwell) + (0.5) * (scalar_v2_cancwell - scalar_v1_cancwell) * ((((2.0/3.0) * (1.0 / (2.0 * omega * (L / (2.0 * R)))) * log((cosh(omega * (((z - midpoint_cancwell) + (L/2.0)) / R)))/(cosh(omega * (((z - midpoint_cancwell) - (L/2.0)) / R))))) + ((1.0/3.0) * (4.0 / (M_PI * M_PI * omega * ((L / 2.0) / R))) * atan((sinh((M_PI * omega * (z - midpoint_cancwell)) / (2.0 * R))) / (cosh((M_PI * omega * (L / 2.0)) / (2.0 * R)))) * atan((sinh((M_PI * omega * (L / 2.0)) / (2.0 * R))) / (cosh((M_PI * omega * (z - midpoint_cancwell)) / (2.0 * R)))) * cosh((M_PI * omega * (z - midpoint_cancwell)) / (2.0 * R)))) + ((((x * x + y * y) * omega) / (2.0 * R * R)) * (1.0 / cosh((omega * (z - midpoint_cancwell)) / R)) * (1.0 / cosh((omega * (z - midpoint_cancwell)) / R)) * tanh((omega * (z - midpoint_cancwell)) / R)));
            // }
            // if (fabs((omega / R) * (z - midpoint_cancwell)) > 1.0)
            // {
            //     if (z - midpoint_cancwell > 0.0)
            //     {
            //         potential_cancwell = (0.5) * (scalar_v1_cancwell + scalar_v2_cancwell) + (0.5) * (scalar_v2_cancwell - scalar_v1_cancwell) * ((((2.0/3.0) * (1.0 / (2.0 * omega * (L / (2.0 * R)))) * log((cosh(omega * (((z - midpoint_cancwell) + (L/2.0)) / R)))/(cosh(omega * (((z - midpoint_cancwell) - (L/2.0)) / R))))) + ((1.0/3.0) * (1.0 - ((4.0 / (M_PI * M_PI * omega * ((L / 2.0) / R))) * atan(1.0 / (sinh((M_PI * omega * (z - midpoint_cancwell)) / (2.0 * R)))) * atan((sinh((M_PI * omega * (L / 2.0)) / (2.0 * R))) / (cosh((M_PI * omega * (z - midpoint_cancwell)) / (2.0 * R)))) * cosh((M_PI * omega * (z - midpoint_cancwell)) / (2.0 * R)))))) + ((((x * x + y * y) * omega) / (2.0 * R * R)) * (1.0 / cosh((omega * (z - midpoint_cancwell)) / R)) * (1.0 / cosh((omega * (z - midpoint_cancwell)) / R)) * tanh((omega * (z - midpoint_cancwell)) / R)));
            //     }
            //     if (z - midpoint_cancwell < 0.0)
            //     {
            //         potential_cancwell = (0.5) * (scalar_v1_cancwell + scalar_v2_cancwell) + (0.5) * (scalar_v2_cancwell + scalar_v1_cancwell) * ((((2.0/3.0) * (1.0 / (2.0 * omega * (L / (2.0 * R)))) * log((cosh(omega * (((-1.0 * (z - midpoint_cancwell)) + (L/2.0)) / R)))/(cosh(omega * (((-1.0 * (z - midpoint_cancwell)) - (L/2.0)) / R))))) + ((1.0/3.0) * (1.0 - ((4.0 / (M_PI * M_PI * omega * ((L / 2.0) / R))) * atan(1.0 / (sinh((M_PI * omega * (-1.0 * (z - midpoint_cancwell))) / (2.0 * R)))) * atan((sinh((M_PI * omega * (L / 2.0)) / (2.0 * R))) / (cosh((M_PI * omega * (-1.0 * (z - midpoint_cancwell))) / (2.0 * R)))) * cosh((M_PI * omega * (-1.0 * (z - midpoint_cancwell))) / (2.0 * R)))))) - ((((x * x + y * y) * omega) / (2.0 * R * R)) * (1.0 / cosh((omega * (-1.0 * (z - midpoint_cancwell))) / R)) * (1.0 / cosh((omega * (-1.0 * (z - midpoint_cancwell))) / R)) * tanh((omega * (-1.0 * (z - midpoint_cancwell))) / R)));
            //     }
            // }

            // Update Counters, Origin //
            potential += (potential_affwell + potential_cancwell);
            midpoint_affwell += (L + d);
            midpoint_cancwell += (L + d);
        }
        return potential;
    }

    // Electric Field Calculation Subroutine // 
    void VectorPotentialCalculation(double x, double y, double z, double *Ex, double *Ey, double *Ez)
    {
        // Define Electric Field Vector Components //
        double deltaVx, deltaVy, deltaVz;
        double dev = 0.001;

        // Run the Scalar Approximation Subroutine for x, y, and z //
        deltaVx = ScalarPotentialApproximation(x + dev, y, z) - ScalarPotentialApproximation(x - dev, y, z);
        deltaVy = ScalarPotentialApproximation(x, y + dev, z) - ScalarPotentialApproximation(x, y - dev, z);
        deltaVz = ScalarPotentialApproximation(x, y, z + dev) - ScalarPotentialApproximation(x, y, z - dev);

        // Calculate E-Field Components //
        *Ex = -deltaVx / (2.0 * dev);
        *Ey = -deltaVy / (2.0 * dev);
        *Ez = -deltaVz / (2.0 * dev);
    }

/////////////////////////////////
// Saved Parameter Subroutines //
/////////////////////////////////

    // Parameter Saving Subroutine //
    void SaveParameters(state *initial, state *proton)
    {
        // Compute Values to be Saved //
        double xi = initial->x, yi = initial->y, zi = initial->z;
        double xf = proton->x, yf = proton->y, zf = proton->z;

        double cyclo_radi = CyclotronRadius(initial->vx, initial->vy, initial->vz);
        double cyclo_radf = CyclotronRadius(proton->vx, proton->vy, proton->vz);

        double azimuthali = AzimuthalAngle(initial->x, initial->y);
        double magtron_phase = proton->magtron_phase;

        double bounces = proton->bounces;

        double Ei = TotalEnergy(initial->x, initial->y, initial->z, initial->vx, initial->vy, initial->vz);
        double Ef = TotalEnergy(proton->x, proton->y, proton->z, proton->vx, proton->vy, proton->vz);

        // Output Values //
        printf("Initial Position: x = %f, y = %f, z = %f\n", xi, yi, zi);
        printf("Final Position: x = %f, y = %f, z = %f\n", xf, yf, zf);
        printf("Initial Energy: E = %f\n", Ei);
        printf("Final Energy: E = %f\n", Ef);
        printf("Total Bounces: %f\n", bounces);
        printf("Cyclotron Radius: initial radius = %f, final radius = %f\n", cyclo_radi, cyclo_radf);
        printf("Magnotron Phase: %f\n", magtron_phase);

        return;
    }

    // Cyclotron Radius Calculation Subroutine // 
    double CyclotronRadius(double vx, double vy, double vz)
    {
        double cyc_r = (mp * sqrt(vx * vx + vy * vy + vz * vz)) / (q * sqrt(Bx * Bx + By * By + Bz * Bz));
        return cyc_r;
    }

    // Compute Azimuthal Angle //
    double AzimuthalAngle(double x, double y)
    {
        double azimuthal_angle = atan2(y, x);
        return azimuthal_angle;
    }

    // Track Magnetron Motion through the Azimuthal Angle //
    void MagnetronMotion(double theta_placeholder, state *proton)
    {
        double theta = AzimuthalAngle(proton->x, proton->y);
        double delta_theta = theta - theta_placeholder;
        
        if (delta_theta > M_PI)
        {
            delta_theta -= 2 * M_PI;
        }
        if (delta_theta < -M_PI)
        {
            delta_theta += 2 * M_PI;
        }

        proton->magtron_phase += delta_theta;
    }

    // Bounce Checker //
    void BounceChecker(double vz_placeholder, state *proton)
    {
        if ((vz_placeholder <= 0 && proton->vz > 0) || (vz_placeholder >= 0 && proton->vz < 0))
        {
            proton->bounces += 1;
        }
    }

////////////////////////////
// Randomizer Subroutines //
////////////////////////////

    // Ranged Randomizer Subroutine //
    double RangedRandomizer(double min, double max) 
    {
        double range = (max - min); 
        double div = RAND_MAX / range;
        return min + (rand() / div);
    }

    // Randomized Position Subroutine //
    void PositionRandomizer(double *xi, double *yi, double *zi, double *ri)
    {
        double TotTrapLength = (16 * L) + (15 * d);

        *ri = sqrt(RangedRandomizer(0, 1)) * 0.04;
        double theta_cyl = RangedRandomizer(0, (M_PI * 2));
        *xi = *ri * cos(theta_cyl);
        *yi = *ri * sin(theta_cyl);
        *zi = RangedRandomizer(0, TotTrapLength);
    }

    // Randomized Velocity Subroutine //
    void VelocityRandomizer(double *vxi, double *vyi, double *vzi, double initial_KE)
    {
        double rho = sqrt((2 * initial_KE * q) / mp);
        double theta_sph = RangedRandomizer(0, (2 * M_PI));
        double cos_phi_sph = RangedRandomizer(-1, 1);

        *vxi = rho * cos(theta_sph) * sin(acos(cos_phi_sph));
        *vyi = rho * sin(theta_sph) * sin(acos(cos_phi_sph));
        *vzi = rho * cos_phi_sph;
    }

///////////////////////////////
// Miscellaneous Subroutines //
///////////////////////////////

    // Total Energy (eV) Calculation Subroutine //
    double TotalEnergy (double x, double y, double z, double vx, double vy, double vz)
    {
        double KE, PE;
        KE = (1.0/2.0) * mp * (vx * vx + vy * vy + vz * vz) / q;
        PE = ScalarPotentialApproximation(x, y, z) * q;

        return KE + PE;
    }