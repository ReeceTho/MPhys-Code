#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dlfcn.h>
#include <sys/wait.h> 
#include "num_in.h"
#include "num_out.h"
#include "VandP.h"
#include "dynamic_cs.h"
#include "rootDir.h" 
#include <time.h>
#define MAX_LINE_LENGTH 1024 

double rnd_slg(double xmin, double xmax) {
    //randomly scans with a log scale as well as a negative side
    double sgn = 1.0;
    if (xmin < 0) {
        xmin = 1E-10;
        sgn = (rand() % 2 == 0) ? 1 : -1; // Randomly return +1 or -1 
    }

    double log_min = log(xmin);
    double log_max = log(xmax);

    // Generate a random value on the log scale
    double random_log = log_min + ((double)rand() / RAND_MAX) * (log_max - log_min);

    return sgn * exp(random_log);
}

double rnd_lin(double xmin, double xmax) {
    return xmin + (xmax - xmin) * ((double)rand() / RAND_MAX);
}

int main(void)
{
    int err, i = 0;

    /* INTPUT PARAMETERS (to scan over) */

    /* OUTPUT PARAMETERS */
    double wh, braa;
    txtList branchings_MP, branchings_HD;

    // Set model dir here
    char mdldir[] = "models";

    // Set model number and number of points to collect, mdlnr is your model number
    int mdlnr = 6;

    // a model to switch between to reset values when reloading
    setModel(mdldir, mdlnr);

    srand(time(NULL)); // seed random number generator

    // Try to remove old output file if exists
    if (remove("scan2.dat") == -1)
        perror("Error in deleting a file");

    // Open the output file for appending
    FILE* file2 = fopen("scan2.dat", "a+"); // Open once before the loop
    if (file2 == NULL) {
        printf("Error: Could not open scan2.dat\n");
        return 1;
    }

    // Writing header to output file
    fprintf(file2, "MD1 \t MD2\t DMP \t DM2 \t DM3 \t Br(h2->W+h-) \t Br(h2->W-h+) \t Br(H->Z h1) \t Br(h2->e+e-) \t Br(h2->mu+mu-) \t Br(h2->tau+tau-) \t Br(h2->neutrino+ neutrino-) \t Br(h2->mu_neutrino+ mu_neutrino-) \t Br(h2-> tau-lepton+ tau_lepton-)\n");

    // Open the scan_tot.dat file for reading
    FILE* file = fopen("scan_tot.dat", "r");
    if (file == NULL) {
        printf("Error: Could not open scan_tot.dat\n");
        fclose(file2); // Close the file2 before exiting
        return 1;
    }

    // Read and ignore the header line
    char line[MAX_LINE_LENGTH];
    fgets(line, sizeof(line), file);

    // Start randomizing loop
    while (fgets(line, sizeof(line), file) != NULL) {
        // Remove newline character at the end of the line
        line[strcspn(line, "\n")] = 0;

        double MD1, DMP, DM2, DM3, MD2;

        // Parse the line for MD1, DMP, and DM3
        if (sscanf(line, "%lf %lf %lf", &MD1, &DMP, &DM3) == 3) {
            DM2 = DM3 + DMP;
            MD2 = DM2 + MD1;

            // Debugging: Print parsed values
            printf("Line %d: MD1=%.6f, DMP=%.6f, DM2=%.6f, DM3=%.6f, MD2=%.6f\n", i, MD1, DMP, DM2, DM3, MD2);


            // Reset model to ensure correct recalculation
            setModel(mdldir, mdlnr);

            // Assigning values to the model parameters
            err = assignValW("MD1", MD1);
            err = assignValW("DMP", DMP);
            err = assignValW("DM3", DM3);

            // Calculate public constraints
            err = calcMainFunc();
            if (err != 0) {
                printf("Can not calculate constrained parameter %s\n", varNames[err]);
                i--;
            }
            else {
                // Collect more output values if the point survives constraints
                double wHD = pWidth("~h2", &branchings_HD);
                double BrHD__Wphm = findBr(branchings_HD, "W+,~h-");
                double BrHD__Wmhp = findBr(branchings_HD, "W-,~h+");
                double BrHD__Zh1 = findBr(branchings_HD, "Z,~h1");
                double BrHD__electron = findBr(branchings_HD, "e1,E1,~h1");
                double BrHD__muon = findBr(branchings_HD, "e2,E2,~h1");
                double BrHD__tau_lepton = findBr(branchings_HD, "e3,E3,~h1");
                double BrHD__neutrino = findBr(branchings_HD, "n1, N1, ~h1");
                double BrHD__mu_neutrino = findBr(branchings_HD, "n2, N2, ~h1");
                double BrHD__tau_neutrino = findBr(branchings_HD, "n3, N3, ~h1");

                // Write the output to file
                fprintf(file2, "%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \t%e \n",
                    MD1, MD2, DMP, DM2, DM3, BrHD__Wphm, BrHD__Wmhp, BrHD__Zh1,
                    BrHD__electron, BrHD__muon, BrHD__tau_lepton,
                    BrHD__neutrino, BrHD__mu_neutrino, BrHD__tau_neutrino);
            }
            i++;
        }
        else {
            printf("Warning: Skipped line %d due to formatting error\n", i);
        }
    }

    // Close files
    fclose(file);
    fclose(file2);
    return 0;
}
