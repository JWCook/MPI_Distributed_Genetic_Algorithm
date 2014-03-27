/* ========================================================================== */
/* Functions used for initializing a starting population. Includes importing  */
/* and exporting parameters to and from files.                                */
/* ========================================================================== */
#include <errno.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "init.h"
#include "fitness.h"
#include "mt_mpi.h"
#include "report.h"
#include "types.h"


/* -------------------------------------------------------------------------- */
/* Initializes a population at generation 0.                                  */
/* Program parameters may be initialized to default values, user-specified    */
/* values, or values from a file. Initial population is randomly generated.   */
/* This function is responsible for initializing the random number generator  */
/* with a seed value (if specified).                                          */
/* subpop        : The poulation and other parameters to initialize           */
/* -------------------------------------------------------------------------- */
void init_population(deme *subpop, int argc, char *argv[]) {
    int i, j, my_rank, n_procs, init_type = 0;
    char *filename = "";
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    
    subpop->fit_tot = 0.0;
    subpop->fit_avg = 0.0;
    subpop->fit_prev = 0.0;
    subpop->fit_novar = 0;
    subpop->fit_max = 0;
    subpop->fit_min = 0;
    subpop->cur_gen = 0;
    subpop->complete = 0;

    // Get any command line arguments
    if        (argc == 1)                                    init_type = 0;
    else if    (argc == 2 && strcmp(argv[1], "-d")  == 0)    init_type = 0;
    else if    (argc == 2 && strcmp(argv[1], "-u")  == 0)    init_type = 1;
    else if    (argc == 3 && strcmp(argv[1], "-de") == 0)    init_type = 2;
    else if    (argc == 3 && strcmp(argv[1], "-ue") == 0)    init_type = 3;
    else if    (argc == 3 && strcmp(argv[1], "-i")  == 0)    init_type = 4;
    else usage();
    if (init_type >= 2 && init_type <= 4)                    filename = argv[2];
    
    // Quick fix; TODO: initialize on processor 0 then send data to rest
    if((init_type == 1 || init_type == 3) && n_procs > 1) {
        if (my_rank == 0)
        printf("Please run on one processor when initializing from input\n"); 
        MPI_Finalize();
        exit(-1);
    }
    
    // Set program parameters (from either defaults, user, or file)
    if (init_type == 1 || init_type == 3) get_input(subpop);
    else if (init_type == 0) {
        subpop->rand_seed     = DEFAULT_RAND_SEED;
        subpop->pop_size     = DEFAULT_POP_SIZE;
        subpop->end_type     = DEFAULT_END_TYPE;
        subpop->end_gen         = DEFAULT_END_GENERATION;
        subpop->ff_type         = DEFAULT_FF_TYPE;
        subpop->conv_gen     = DEFAULT_CONV_GENS;
        subpop->conv_var     = DEFAULT_CONV_VARIATION;
        if (subpop->ff_type == FF_SIMPLE)
            subpop->f_thresh = DEFAULT_F_THRESH_SIMPLE;
        else if (subpop->ff_type == FF_SHPATH)
            subpop->f_thresh = DEFAULT_F_THRESH_SHPATH;
        else subpop->f_thresh = 0;
    }
    if (init_type == 4) import_population(subpop, filename);

    // Allocate and randomize population
    mt_init(my_rank, subpop->rand_seed);
    if         (subpop->ff_type == FF_SIMPLE)   subpop->chr_size = CHR_SIZE_SIMPLE;
    else if (subpop->ff_type == FF_SHPATH)      subpop->chr_size = CHR_SIZE_SHPATH;
    else                                        subpop->chr_size = 0;
    subpop->old_pop = (org**) malloc(sizeof(org*) * subpop->pop_size);
    subpop->new_pop = (org**) malloc(sizeof(org*) * subpop->pop_size);
    for (i = 0; i < subpop->pop_size; i++) {
        subpop->old_pop[i] = (org*) malloc(sizeof(org));
        subpop->old_pop[i]->chr = (char*) malloc(sizeof(char)
            * subpop->chr_size);
        subpop->new_pop[i] = (org*) malloc(sizeof(org));
        subpop->new_pop[i]->chr = (char*) malloc(sizeof(char) 
            * subpop->chr_size);
        for (j = 0; j < subpop->chr_size; j++)
            subpop->old_pop[i]->chr[j] = mt_rand_bit(my_rank);
    }
    
    if(init_type == 1 || init_type == 3 || init_type == 4) test_input(subpop);
    if(init_type == 2 || init_type == 3) export_population(subpop, filename);
    else fitness(subpop);
}


/* -------------------------------------------------------------------------- */
/* Imports program parameters from a file following the format specfied in    */
/* export_population.                                                         */
/* -------------------------------------------------------------------------- */
void import_population(deme *subpop, char *filename) {
    int i;
    char varname1[BUFFER_SIZE];
    char varname2[BUFFER_SIZE];
    FILE *fp;

    errno = 0;
    fp = fopen(filename, "r");
    if (errno != 0) {
        fprintf(stderr, "Error %i: %s\n", errno, strerror(errno));
        exit(-1);
    }
    
    subpop->rand_seed     = get_value(fp, "rand_seed:");
    subpop->pop_size     = get_value(fp, "pop_size:");
    subpop->ff_type         = get_value(fp, "ff_type:");
    subpop->end_type     = get_value(fp, "end_type:");
    if (subpop->end_type == 0)
        subpop->end_gen     = get_value(fp, "end_gen:");
    else if (subpop->end_type == 1 || subpop->end_type == 2)
        subpop->f_thresh = get_value(fp, "f_thresh:");
    else if (subpop->end_type == 3) {
        subpop->conv_gen = get_value(fp, "conv_gen:");
        subpop->conv_var = get_value(fp, "conv_var:");
    }
    
    if (subpop->ff_type == FF_SHPATH) {
        subpop->s          = get_point(fp, "s:");
        subpop->t         = get_point(fp, "t:");
        subpop->n_objs     = get_value(fp, "n_objs:");
        subpop->objs = (object**) malloc(sizeof(object*) * subpop->n_objs);
        for (i = 0; i < subpop->n_objs; i++) {
            sprintf(varname1, "obj_%i_center:", i);
            sprintf(varname2, "obj_%i_radius:", i);
            subpop->objs[i] = (object*) malloc(sizeof(object));
            subpop->objs[i]->center = get_point(fp, varname1);
            subpop->objs[i]->radius = get_value(fp, varname2);
        }
    }
    
    fclose(fp);
}


/* -------------------------------------------------------------------------- */
/* Exports program parameters to a file. Format:                              */
/*        #Header                                                             */
/*        rand_seed: <value>                                                  */
/*        pop_size: <value>                                                   */
/*        ff_type: <value>                                                    */
/*        end_type: <value>                                                   */
/*        *end_gen: <value>                                                   */
/*        *f_thresh: <value>                                                  */
/*        *conv_gen: <value>                                                  */
/*        *conv_var: <value>                                                  */
/*        *s: <values>                                                        */
/*        *t: <values>                                                        */
/*        *n_objs: <value>                                                    */
/*        *obj_0_center: <values>                                             */
/*        *obj_0_radius: <value>                                              */
/*        *...                                                                */
/*        *obj_n_center: <values>                                             */
/*        *obj_n_radius: <value>                                              */
/*                                                                            */
/*        * = optional, depending on other settings                           */
/* -------------------------------------------------------------------------- */
void export_population(deme *subpop, char *filename) {
    int i;
    FILE *fp;
    
    errno = 0;
    fp = fopen(filename, "w");
    if (errno != 0) {
        fprintf(stderr, "Error %i: %s\n", errno, strerror(errno));
        exit(-1);
    }

    fprintf(fp, "#File %s generated by ga/init\n", filename);
    fprintf(fp, "rand_seed: %i\n", subpop->rand_seed);
    fprintf(fp, "pop_size: %i\n", subpop->pop_size);
    fprintf(fp, "ff_type: %i\n", subpop->ff_type);
    fprintf(fp, "end_type: %i\n", subpop->end_type);
    
    if(subpop->end_type == M_FIXED_GENERATIONS)
        fprintf(fp, "end_gen: %i\n", subpop->end_gen);
    else if(subpop->end_type == M_MAX_FITNESS_THRESHHOLD
         || subpop->end_type == M_AVG_FITNESS_THRESHHOLD)
        fprintf(fp, "f_thresh: %.1f\n", subpop->f_thresh);
    else if(subpop->end_type == M_LOCAL_CONVERGENCE) {
        fprintf(fp, "conv_gen: %i\n", subpop->conv_gen);
        fprintf(fp, "conv_var: %.1f\n", subpop->conv_var);
    }
    
    if (subpop->ff_type == FF_SHPATH) {
        fprintf(fp, "s: (%i,%i,%i)\n", subpop->s->x,subpop->s->y,subpop->s->z);
        fprintf(fp, "t: (%i,%i,%i)\n", subpop->t->x,subpop->t->y,subpop->t->z);
        fprintf(fp, "\nn_objs: %i\n", subpop->n_objs);
        for (i = 0; i < subpop->n_objs; i++) {
            fprintf(fp, "obj_%i_center: (%i,%i,%i)\n", i,
                    subpop->objs[i]->center->x,
                    subpop->objs[i]->center->y,
                    subpop->objs[i]->center->z);
            fprintf(fp, "obj_%i_radius: %i\n", i, subpop->objs[i]->radius);
        }
    }
    
    printf("File %s successfully written\n", filename);
    fclose(fp);
    exit(1);
}


/* -------------------------------------------------------------------------- */
/* Collects user input for program parameters                                 */
/* -------------------------------------------------------------------------- */
void get_input(deme *subpop) {
    int i;

    printf("Enter population size: ");
    subpop->pop_size = get_value(stdin, NULL);
    
    printf("Enter fitness function type (0-1): ");
    subpop->ff_type = get_value(stdin, NULL);
    
    printf("Enter termination type (0-3): ");
    subpop->end_type = get_value(stdin, NULL);
    
    if (subpop->end_type == 0) {
        printf("Enter end generation: ");
        subpop->end_gen = get_value(stdin, NULL);
    }
    else if (subpop->end_type == 1 || subpop->end_type == 2) {
        printf("Enter fitness threshhold: ");
        subpop->f_thresh = get_value(stdin, NULL);
    }
    else if (subpop->end_type == 3) {
        printf("Enter convergence generation threshhold: ");
        subpop->conv_gen = get_value(stdin, NULL);

        printf("Enter convergence variation: ");
        subpop->conv_var = get_value(stdin, NULL);
    }
    
    printf("Enter RNG seed (0 to generate one): ");
    subpop->rand_seed = get_value(stdin, NULL);
    
    if (subpop->ff_type == FF_SHPATH) {
        printf("Enter starting point in the format (x,y,z): ");
        subpop->s = get_point(stdin, NULL);
        
        printf("Enter end point in the format (x,y,z): ");
        subpop->t = get_point(stdin, NULL);
        
        printf("Enter number of objects ");
        subpop->n_objs = get_value(stdin, NULL);
        subpop->objs = (object**) malloc(sizeof(object*) * subpop->n_objs);
        
        for (i = 0; i < subpop->n_objs; i++) {
            subpop->objs[i] = (object*) malloc(sizeof(object));
            printf("Enter object %i center in the format (x,y,z): ", i);
            subpop->objs[i]->center = get_point(stdin, NULL);
            
            printf("Enter object %i radius: ", i);
            subpop->objs[i]->radius = get_value(stdin, NULL);
        }
    }
}


/* -------------------------------------------------------------------------- */
/* Validate input received from user or read from input file                  */
/* -------------------------------------------------------------------------- */
void test_input(deme* subpop) {
    int i, is_invalid = 0;
    
    if (subpop->pop_size < 1) {
        fprintf(stderr, "Error: Invalid value for population\n");
        is_invalid = 1;
    }
    else if (subpop->pop_size > 10000) 
        fprintf(stderr, "Warning: excessively large population\n");
        
    if (subpop->ff_type < 0 || subpop->ff_type > 1) {
        fprintf(stderr, "Error: Invalid value for fitness function type\n");
        is_invalid = 1;
    }
    
    if (subpop->end_type < 0 || subpop->end_type > 3) {
        fprintf(stderr, "Error: Invalid value for termination type\n");
        is_invalid = 1;
    }

    if (subpop->end_type == 0 && subpop->end_gen < 1) {
        fprintf(stderr, "Error: Invalid value for end generation\n");
        is_invalid = 1;
    }
    else if ((subpop->end_type == 1 || subpop->end_type == 2)
            && subpop->f_thresh < 1) {
        fprintf(stderr, "Error: Invalid value for fitness threshhold\n");
        is_invalid = 1;
    }
    else if (subpop->end_type == 3 
            && (subpop->conv_gen < 1 || subpop->conv_var < 1)) {
        fprintf(stderr, "Error: Invalid convergence values\n");
        is_invalid = 1;
    }
    
    if (subpop->ff_type == FF_SHPATH) {
        if (!valid_loc(subpop->s)) {
            fprintf(stderr, "Error: Invalid start point\n");
            is_invalid = 1;
        }
        
        if (!valid_loc(subpop->t)) {
            fprintf(stderr, "Error: Invalid start point\n");
            is_invalid = 1;
        }
    
        for (i = 0; i < subpop->n_objs; i++) {
            if (subpop->objs[i]->center == NULL) {
                fprintf(stderr, "Error: No entry for object %i\n",i);
                is_invalid = 1;
            }
            else if (!valid_loc(subpop->objs[i]->center)) {
                fprintf(stderr,"Error: Invalid coordinates for object %i\n",i);
                is_invalid = 1;
            }
            
            if (subpop->objs[i]->center!=NULL && subpop->objs[i]->radius < 1) {
                fprintf(stderr, "Error: Invalid radius for object %i\n",i);
                is_invalid = 1;
            }
            
            if (pt_dist(subpop->s, subpop->objs[i]->center) 
                < subpop->objs[i]->radius) {
                fprintf(stderr,"Error: object %i overlaps start point\n",i);
                is_invalid = 1;
            }
            
            if (pt_dist(subpop->t, subpop->objs[i]->center) 
                < subpop->objs[i]->radius) {
                fprintf(stderr,"Error: object %i overlaps end point\n",i);
                is_invalid = 1;
            }
        }
    }
    
    if (is_invalid) usage2();
}


/* -------------------------------------------------------------------------- */
/* Get a value from the specified file as a floating point number             */
/* fp            : A pointer to the stream being read                         */
/* varname        : The name of the value to read, if reading from a file     */
/* return        : The number value immediately following varname in fp, if   */
/*                  any; 0 otherwise                                          */
/* -------------------------------------------------------------------------- */
double get_value(FILE *fp, char *varname) {
    int status;
    char token[BUFFER_SIZE];
    
    if (varname != NULL) {
        rewind(fp);
        do status = fscanf(fp, "%s", token);
        while (!feof(fp) && strcmp(token, varname) != 0);
    }
    
    if (!feof(fp)) {
        status = fscanf(fp, "%s", token);
        return atof(token);
    }
    else return -1;
}


/* -------------------------------------------------------------------------- */
/* Make a new point structure from a string in the format "(x,y,z)".          */
/* If a string is incorrectly formatted, 0s are used as coordinates           */
/* fp            : A pointer to the stream being read                         */
/* s            : The label of the coordinates to look for                    */
/* -------------------------------------------------------------------------- */
point *get_point(FILE *fp, char *varname) {
    int status;
    char str[BUFFER_SIZE];
    char *x_str, *y_str, *z_str;
    point *p = (point*) malloc(sizeof(point));
    
    if (varname != NULL) {
        rewind(fp);
        do status = fscanf(fp, "%s", str);
        while (!feof(fp) && strcmp(str, varname) != 0);
    }
    
    if        (!feof(fp)) status = fscanf(fp, "%s", str);
    else    return NULL;
        
    x_str = strtok(str,  "( ,)");
    y_str = strtok(NULL, "( ,)");
    z_str = strtok(NULL, "( ,)");
    
    if (x_str != NULL) p->x = atoi(x_str);     else p->x = SP_BOUND+1;
    if (y_str != NULL) p->y = atoi(y_str);     else p->y = SP_BOUND+1;
    if (z_str != NULL) p->z = atoi(z_str);     else p->z = SP_BOUND+1;
    
    return p;
}

