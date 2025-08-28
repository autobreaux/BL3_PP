#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// Global Constants //
    const double Bz = 0.5; // Magnetic field z component //
    const double q = 1.602e-19, mp = 1.673e-27; // Proton charge, mass //
    const double R = 0.035, d = 0.0048; // Cylinder radius, midpoint between cylinders //
    const double L = 0.0324; // Length of Cylinder //
    const double omega = 1.3262275; // Bertam constant //

int main(){ 
    double z = 0.0;
    double x = 0.0, y = 0.0;


    // Open Binary File //
    FILE *data_file = fopen("approx_test.bin", "wb");
    if (data_file == NULL) {
        perror("Error opening file");
        exit(1);
    }

    // Testing Approx. Main Loop //
    while(z <= 0.5904)
    {
        // Cylinder Voltages in Order from Left to Right //
        double voltages[] = {850.0, 850.0, 850.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 850.0, 850.0, 850.0};
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

            // // // Second Consider the Cancelling Well //
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
    

    // Write Into Igor File //
    double data[4] = {z, potential, potential_affwell, potential_cancwell};
    fwrite(data, sizeof(double), 4, data_file);

    z += 0.0005904;
    }
    // Close Igor File //
    fclose(data_file);
}