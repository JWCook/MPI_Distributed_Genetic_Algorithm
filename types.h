#ifndef TYPES_H_
#define TYPES_H_


/* -------------------------------------------------------------------------- */
/* A structure representing a three-dimensional point                         */
/* -------------------------------------------------------------------------- */
typedef struct {
    int x;
    int y;
    int z;
} point;


/* -------------------------------------------------------------------------- */
/* A structure representing a 3-dimensional object, approximated by a sphere  */
/* center        : The center of the sphere                                   */
/* radius        : The radius of the sphere                                   */
/* -------------------------------------------------------------------------- */
typedef struct {
    point    *center;
    int        radius;
} object;


/* -------------------------------------------------------------------------- */
/* A struct representing an individual organism.                              */
/* chr        : The individual's chromosome in the form of a bit string       */
/* fitness    : The individual's estimated fitness level                      */
/* parent1    : The index of one of the individual's parents                  */
/* parent2    : The index of one of the individual's parents                  */
/* -------------------------------------------------------------------------- */
typedef struct {
    double      fitness;
    int         parent1;
    int         parent2;
    char        *chr;
} org;


/* -------------------------------------------------------------------------- */
/* A struct representing a (sub)population.                                   */
/* old_pop      : The members of the current generation                       */
/* new_pop      : The members of the new generation being generated           */
/* rand_seed    : The seed used to initialize the random number generator     */
/* chr_size     : The size of the chromosomes in this population              */
/* fit_tot      : The total fitness of this population                        */
/* fit_avg      : The average fitness of this population                      */
/* fit_max      : The index of the most fit member of this population         */
/* fit_min      : The index of the least fit member of this population        */
/* fit_prev     : The average fitness of the previous generation              */
/* fit_novar    : The number consecutive generations for which this           */
/*                  population's fit_avg has varied by less than conv_var     */
/* cur_gen      : The current generation                                      */
/* pop_size     : The size of this population                                 */
/* ff_type      : The fitness function to use                                 */
/* end_type     : The type of termination condition                           */
/* end_gen      : Number of generations, if end type is M_FIXED_GENERATIONS   */
/* f_thresh     : Average or max fitness threshhold, if end type is either    */
/*                  M_AVG_FITNESS_THRESHHOLD or M_MAX_FITNESS_THRESHHOLD      */
/* conv_gen     : Number of generations to test for a local convergence       */
/* conv_var     : Threshhold for testing if a local convergence has occured,  */
/*                  if end type is M_LOCAL_CONVERGENCE; if the average fitness*/
/*                  varies by less than conv_var for conv_gens generations,   */
/*                  the algorithm has converged on a local min/max            */
/* complete     : A non-zero value flags this population as terminated        */
/*                                                                            */
/* If the shortest path fitness function is being used:                       */
/* n_objs       : The number of obstacles                                     */
/* objs         : The set of obstacles                                        */
/* s            : Start point                                                 */
/* t            : End point                                                   */
/* -------------------------------------------------------------------------- */
typedef struct {
    org         **old_pop;
    org         **new_pop;
    int         rand_seed;
    int         chr_size;
    double      fit_tot;
    double      fit_avg;
    int         fit_max;
    int         fit_min;
    int         fit_prev;
    int         fit_novar;
    int         cur_gen;
    int         pop_size;
    int         ff_type;
    int         end_type;
    int         end_gen;
    double      f_thresh;
    int         conv_gen;
    double      conv_var;
    int         complete;
    
    int         n_objs;
    object      **objs;
    point       *s;
    point       *t;
} deme;


#endif

