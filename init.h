#ifndef INIT_H_
#define INIT_H_
#include "types.h"

void    init_population(deme*, int, char**);
void    import_population(deme*, char*);
void    export_population(deme*, char*);
void    get_input(deme*);
void    test_input(deme*);
double  get_value(FILE*, char*);
point   *get_point(FILE*, char*);

#endif

