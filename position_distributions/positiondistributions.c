#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

// Global Constants //
    const double Bx = 0.0, By = 0.0, Bz = 0.5; // Magnetic field components //
    const double q = 1.602e-19, mp = 1.673e-27; // Proton charge, mass //
    const double R = 0.035, d = 0.0048; // Cylinder radius, distance between cylinders //
    const double L = 0.0324; // Length of cylinder //
    const double omega = 1.3262275; // Bertam constant //
    const double t0 = 0; // Initial time //
    const double Ex0 = 0, Ey0 = 0, Ez0 = 0; // Initial Electric Field Components //
    double z1_maxpotential = 0.0371952, z2_maxpotential = 0.520733; // Z-coords of max potentials //

// Ranged Randomizer Subroutine //
double RangedRandomizer(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Randomized Position Subroutine //
void PositionRandomizer(double *x0, double *y0, double *z0, double *r0)
{
    double TotTrapLength = (16 * L) + (15 * d);

    *r0 = sqrt(RangedRandomizer(0, 1)) * 0.04;
    double theta_cyl = RangedRandomizer(0, (M_PI * 2));
    *x0 = *r0 * cos(theta_cyl);
    *y0 = *r0 * sin(theta_cyl);
    *z0 = RangedRandomizer(0, TotTrapLength);

}

// Main Function //
int main()
{
    // Runtime Declarations, Constants //
    double x0, y0, z0, r0;
    double vx0, vy0, vz0;

    // Open Binary File(s) //
    FILE *file1 = fopen("randompositions.bin", "wb");
    if (file1 == NULL) 
    {
        perror("Error opening file");
        exit(1);
    }

    // Begin Timer Randomizer //
    srand(time(NULL));

    // Main Loop //
    for(int i = 0; i < 10000; i++)
    {
        // Randomized Position Initialization //
        PositionRandomizer(&x0, &y0, &z0, &r0);

        // Write into Binary File //
        double data[4] = {x0, y0, z0, r0};
        fwrite(data, sizeof(double), 4, file1);
    }
    
    // Close Igor File //
    fclose(file1);

    return 0;
}