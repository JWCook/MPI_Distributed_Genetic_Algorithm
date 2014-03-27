#ifndef FITNESS_H_
#define FITNESS_H_
#include "types.h"

void 	fitness(deme*);
void	fitness_simple(deme*);
void	fitness_shpath(deme*);
double	pt_dist(point*, point*);
int		collision(point*, point*, object*);
int		valid_loc(point*);
void	pt_copy(point*, point*);
point   **make_path(char*, point*, point*);
void    free_path(point**);
int 	binToDecimal(char*, unsigned int, unsigned int);

#endif

