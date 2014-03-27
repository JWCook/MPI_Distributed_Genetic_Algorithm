#ifndef GA_H_
#define GA_H_
#include "types.h"

void	migration(deme*);
int 	selection(deme*);
void 	reproduction(deme*);
void 	crossover(deme*);
void 	mutation(deme*);
void 	check_complete(deme*);
void	sync_complete(deme*);

#endif

