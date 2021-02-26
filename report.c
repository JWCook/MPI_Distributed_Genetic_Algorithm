/* ========================================================================== */
/* Functions for reporting population statistics and other info               */
/* ========================================================================== */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "config.h"
#include "fitness.h"
#include "types.h"
#include "report.h"


/* -------------------------------------------------------------------------- */
/* Print data on an individual member in a readable format                    */
/* -------------------------------------------------------------------------- */
void report_member(deme *subpop, char *chr) {
    int i;

    if (subpop->ff_type == FF_SIMPLE) {
        printf("[x=%03i, y=%03i, a=%03i, b=%03i]",
            binToDecimal(chr, 0,  7),
            binToDecimal(chr, 8,  15),
            binToDecimal(chr, 16, 23),
            binToDecimal(chr, 24, 31));
    }
    else {
        point **path = make_path(chr, subpop->s, subpop->t);
        for (i = 0; i < N_POINTS+2; i++)
            printf("%i: (%i,%i,%i)\n", i, path[i]->x, path[i]->y, path[i]->z);
        free_path(path);
    }
}



/* -------------------------------------------------------------------------- */
/* Print overall population data                                              */
/* -------------------------------------------------------------------------- */
void report_all(deme* subpop) {
    int i, my_rank;
    char chr[subpop->chr_size+1];
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (PRNT_RATE > 0 && (subpop->cur_gen <= 1
        || subpop->cur_gen % PRNT_RATE == PRNT_RATE - 1 || subpop->complete)) {
        if (PRNT_DATA && (PRNT_CHRS || PRNT_INFO || PRNT_FITS)) {
            for (i = 0; i < subpop->pop_size; i++) {
                printf("%03i: ", i+1);
                if (PRNT_CHRS) {
                    strncpy(chr, subpop->new_pop[i]->chr, subpop->chr_size);
                    chr[subpop->chr_size] = '\0';
                    printf("%s ", chr);
                }
                if (PRNT_INFO)
                    report_member(subpop, subpop->new_pop[i]->chr);
                if (PRNT_FITS)
                    printf("%.0f", subpop->new_pop[i]->fitness);
                printf("\n");
            }
        }
        if (PRNT_STAT) {
            printf("[Deme %03i][Gen %06i] ", my_rank, subpop->cur_gen);
            printf("Total:%08.0f ", subpop->fit_tot);
            printf("Avg:%06.0f ", subpop->fit_avg);
            printf("Max[%03i]: ", subpop->fit_max+1);
            printf("%.0f\n", subpop->new_pop[subpop->fit_max]->fitness);
        }
    }
}


/* -------------------------------------------------------------------------- */
/* Finds and prints the most fit member across all sub-populations            */
/* -------------------------------------------------------------------------- */
void report_fittest(deme *subpop) {
    int          source, my_rank, n_procs, global_max = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Status    status;
    char        local_fit[subpop->chr_size];
    char        fittest[n_procs][subpop->chr_size+1];
    double      maxima[n_procs];

    // Every process except 0 sends its most fit member to process 0
    if (my_rank != 0) {
        MPI_Send(subpop->new_pop[subpop->fit_max]->chr, subpop->chr_size,
            MPI_CHAR, 0, 50, MPI_COMM_WORLD);
        MPI_Send(&subpop->new_pop[subpop->fit_max]->fitness, 1, MPI_DOUBLE, 0,
         50, MPI_COMM_WORLD);
    }

    // Process 0 collects the most fit member of it and every other process
    else {
        strncpy(fittest[0], subpop->new_pop[subpop->fit_max]->chr,
            subpop->chr_size);
        fittest[0][subpop->chr_size] = '\0';
        maxima[0] = subpop->new_pop[subpop->fit_max]->fitness;

        for (source = 1; source < n_procs; source++) {
            MPI_Recv(local_fit, subpop->chr_size, MPI_CHAR, source,
                50, MPI_COMM_WORLD, &status);
            MPI_Recv(&maxima[source], 1, MPI_DOUBLE, source, 50,
                MPI_COMM_WORLD,
                &status);
            strncpy(fittest[source], local_fit, subpop->chr_size);
            fittest[source][subpop->chr_size] = '\0';
        }

        usleep(50);

        for(source = 0; source < n_procs; source++)
            if (maxima[source] > maxima[global_max]) global_max = source;

        printf("Best solution found:\n");
        if (PRNT_CHRS) printf("%s\n", fittest[global_max]);
        if (PRNT_INFO) report_member(subpop, fittest[global_max]);
        if (PRNT_FITS) printf("\nFitness: %.0f",maxima[global_max]);
        printf("\n");
    }
}


/* -------------------------------------------------------------------------- */
/* Print a usage statement and exit the program                               */
/* -------------------------------------------------------------------------- */
void usage() {
    printf("\n\nUsage: ga [flag] [FILE]\n");
    printf("\t-d: Initialize data with default values\n");
    printf("\t-u: Initialize data with user-specified values\n");
    printf("\t-de: Initialize data with default values; exports data\n");
    printf("\t    to a file without running the program.\n");
    printf("\t-ue: Initialize with user-specified values; exports data\n");
    printf("\t    to a file without running the program.\n");
    printf("\t-i: Initialize data from a previously generated file.\n");
    printf("\tA valid filename must be specified if importing or exporting.");
    printf("\n\n");

    exit(-1);
}


/* -------------------------------------------------------------------------- */
/* Print info on usage of program parameters and exit the program             */
/* -------------------------------------------------------------------------- */
void usage2() {
    printf("\n\nUser-specified values:\n\n");
    printf("Population size:\nA positive integer value for the size ");
    printf("of the population. For parallel execution, this determines ");
    printf("the size of each sub-population.\n\n");
    printf("Fitness function:\nThe type of fitness function to use. ");
    printf("Possible values:\n");
    printf("0: Maximizes a simple function\n");
    printf("1: Finds the shortest 3D path between 2 points that\n");
    printf("   does not collide with a set of objects\n\n");
    printf("Termination type:\nThe condition on which the program will ");
    printf("terminate. Possible values:\n");
    printf("0: Fixed number of generations\n");
    printf("1: Maximum fitness threshhold\n");
    printf("2: Average fitness threshhold\n");
    printf("3: Local convergence\n\n");
    printf("End generation:\nA positive integer value for the number of ");
    printf("generations, if a termination type of fixed generations was ");
    printf("selected\n\n");
    printf("Fitness threshhold:\nA positive value for an accepted fitness ");
    printf("level, if a termination type of average or maximum fitness ");
    printf("threshhold was selected\n\n");
    printf("Convergence generation threshhold:\nA positive integer value ");
    printf("for the number of generations with low variation that should ");
    printf("be considered a local convergence, if a termination type of ");
    printf("local convergence was selected\n\n");
    printf("Convergence variation:\nA positive value for the variation ");
    printf("allowed for a local convergence, if a termination type of local ");
    printf("convergence was selected\n\n");
    printf("RNG seed:\nA value used to initialize the random rumber ");
    printf("generator. 0 may be used to generate a new seed value. Mutliple ");
    printf("runs of the program with the same seed value will yield the ");
    printf("same results.\n\n");
    printf("Start and end points:\nThese define the end points between ");
    printf("which a 3-dimensional path is contructed, entered in the format ");
    printf("(x,y,z). coordinates must be within the bounds ");
    printf("(%i,%i)\n\n", SP_BOUND*-1, SP_BOUND);
    printf("Number of obstacles:\nThese are the obstacles around which a ");
    printf("path must be constructed. Each object is aproximated by a ");
    printf("sphere, with a center entered in the format (x,y,z) and a ");
    printf("positive integer value for the radius.\n");

    exit(-1);
}

