#define PRNT_RATE                   500     // Rate at which to report data
#define PRNT_DATA                   1       // Report pop data while running
#define PRNT_STAT                   1       // Report pop stats
#define PRNT_CHRS                   0       // Report raw chromosome data
#define PRNT_INFO                   1       // Report human-readable data
#define PRNT_FITS                   1       // Report individual fitness
#define BUFFER_SIZE                 64

#define CROSSOVER_RATE              0.8
#define MUTATION_RATE               0.05

#define FF_SIMPLE                   0
#define FF_SHPATH                   1
#define M_FIXED_GENERATIONS         0
#define M_MAX_FITNESS_THRESHHOLD    1
#define M_AVG_FITNESS_THRESHHOLD    2
#define M_LOCAL_CONVERGENCE         3

#define DEFAULT_RAND_SEED           42
#define DEFAULT_FF_TYPE             FF_SIMPLE
#define DEFAULT_END_TYPE            M_AVG_FITNESS_THRESHHOLD
#define DEFAULT_POP_SIZE            320
#define DEFAULT_END_GENERATION      6000
#define DEFAULT_F_THRESH_SIMPLE     125000.0
#define DEFAULT_F_THRESH_SHPATH     10000.0
#define DEFAULT_CONV_GENS           20
#define DEFAULT_CONV_VARIATION      100

#define CHR_SIZE_SIMPLE             32

#define SP_BOUND                    2048
#define N_POINTS                    32
#define COORD_SIZE                  12
#define CHR_SIZE_SHPATH             COORD_SIZE*3*N_POINTS
#define COLLISION_COST              100

