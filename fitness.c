/* ========================================================================= */
/* Fitness functions and associated helper methods                           */
/* ========================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "fitness.h"
#include "types.h"


/* ------------------------------------------------------------------------- */
/* Call the appropriate fitness function                                     */
/* ------------------------------------------------------------------------- */
void fitness(deme *subpop) {
    if      (subpop->ff_type == FF_SIMPLE)  fitness_simple(subpop);
    else if (subpop->ff_type == FF_SHPATH)  fitness_shpath(subpop);
}


/* ------------------------------------------------------------------------- */
/* Simple fitness function that maximizes the equation:                      */
/* f(x,y) = y^2 + x^2 + a^2 + b^2 - 2xy - 2ab                                */
/* (maximum possible fitness of 130050 with 8-bit values)                    */
/* ------------------------------------------------------------------------- */
void fitness_simple(deme *subpop) {
    int i, x, y, a, b;
    double fit;
    org** pop;
    if (subpop->fit_tot < 1) pop = subpop->old_pop;
    else                     pop = subpop->new_pop;

    subpop->fit_prev = subpop->fit_avg;
    subpop->fit_avg = 0.0;
    subpop->fit_tot = 0.0;

    for (i = 0; i < subpop->pop_size; i++) {
        x = binToDecimal(pop[i]->chr, 0, 7);
        y = binToDecimal(pop[i]->chr, 8, 15);
        a = binToDecimal(pop[i]->chr, 16, 23);
        b = binToDecimal(pop[i]->chr, 24, 31);
        fit = (y*y) + (x*x) + (a*a) + (b*b) - (2*x*y) - (2*a*b);
        if (fit < 0) fit = 0;
        pop[i]->fitness = fit;
        subpop->fit_tot += fit;
        if (fit > pop[subpop->fit_max]->fitness) subpop->fit_max = i;
        if (fit < pop[subpop->fit_min]->fitness) subpop->fit_min = i;

    }

    subpop->fit_avg = subpop->fit_tot / subpop->pop_size;
}


/* ------------------------------------------------------------------------- */
/* Evaluates each population member as a 3D path.                            */
/* Fitness is determined by path length and number of object collisions.     */
/* ------------------------------------------------------------------------- */
void fitness_shpath(deme *subpop) {
    int i, j, k, colls;
    double cost, fit, dist;
    double max_cost = SP_BOUND*40;

    org** pop;
    if (subpop->fit_tot < 1) pop = subpop->old_pop;
    else                     pop = subpop->new_pop;
    point **path;

    subpop->fit_prev = subpop->fit_avg;
    subpop->fit_avg = 0.0;
    subpop->fit_tot = 0.0;

    for (i = 0; i < subpop->pop_size; i++) {
        // Construct and find length of path
        path = make_path(pop[i]->chr, subpop->s, subpop->t);
        dist = 0;
        for (j = 0; j < N_POINTS+1; j++)
            dist += pt_dist(path[j], path[j+1]);

        // Detect collisions between path and obstacles
        colls = 0;
        for (j = 0; j < N_POINTS+1; j++) {
            for (k = 0; k < subpop->n_objs; k++) {
                if (collision(path[j], path[j+1], subpop->objs[k])) colls++;
            }
        }

        // Calculate fitness
        cost = dist + (colls*COLLISION_COST);
        fit = max_cost - cost;
        if (fit < 0) fit = 1;
        pop[i]->fitness = fit;
        subpop->fit_tot += fit;
        if (fit > pop[subpop->fit_max]->fitness) subpop->fit_max = i;
        if (fit < pop[subpop->fit_min]->fitness) subpop->fit_min = i;

        free_path(path);
    }

    subpop->fit_avg = subpop->fit_tot / subpop->pop_size;
    //printf("FIT DONE\n");
}


/* ------------------------------------------------------------------------- */
/* Find the distance between two points                                      */
/* ------------------------------------------------------------------------- */
double pt_dist(point *p1, point *p2) {
    double x = (p1->x - p2->x);
    double y = (p1->y - p2->y);
    double z = (p1->z - p2->z);
    return sqrt((x*x) + (y*y) + (z*z));
}


/* ------------------------------------------------------------------------- */
/* Tests if a line collides with an object.                                  */
/* Since the object is approximated by a sphere, take the perpendicular      */
/* distance from the center of the sphere to the line. If this is within the */
/* object's radius, the line intersects the object.                          */
/* p1           : The start point of the line to test                        */
/* p2           : The end point of the line to test                          */
/* obj          : The object to test                                         */
/* ------------------------------------------------------------------------- */
int collision(point *A, point *B, object *obj) {
    int i1, j1, k1, i2, j2, k2, i3, j3, k3;
    double area, len, dist;

    // Set up vectors AB and AP, where P is the object's center
    i1 = B->x - A->x;           i2 = obj->center->x - A->x;
    j1 = B->y - A->y;           j2 = obj->center->y - A->y;
    k1 = B->y - A->y;           k2 = obj->center->z - A->z;

    // Find cross product AB X AP
    i3 = (j1*k2) - (k1*j2);
    j3 = (k1*i2) - (i1*k2);
    k3 = (i1*j2) - (j1*i2);

    // Find distance from AB to P by (area of AB X AP)/(length of AB)
    area = sqrt((i3*i3) + (j3*j3) + (k3*k3));
    len  = sqrt((i1*i1) + (j1*j1) + (k1*k1));
    dist = area / len;

    // If this is within the object's radius, AB collides with obj
    return (dist < obj->radius);
}


/* ------------------------------------------------------------------------- */
/* Test if a point is in a valid location; e.g., within the defined bounds   */
/* ------------------------------------------------------------------------- */
int valid_loc(point *p) {
    if (p->x > SP_BOUND || p->x < (0-SP_BOUND)|| p->y > SP_BOUND
         || p->y < (0-SP_BOUND)|| p->z > SP_BOUND || p->z < (0-SP_BOUND))
         return 0;
    else return 1;
}


/* ------------------------------------------------------------------------- */
/* Copy the coordinates of one point into another                            */
/* ------------------------------------------------------------------------- */
void pt_copy(point *dest, point *src) {
    dest->x = src->x;
    dest->y = src->y;
    dest->z = src->z;
}


/* ------------------------------------------------------------------------- */
/* Construct a path from a set of points represented by a binary string      */
/* chr          : A binary string representing a list of points. Each        */
/*                coordinate is COOR_SIZE bits long.                         */
/* s            : The start point of the path                                */
/* t            : The end point of the path                                  */
/* ------------------------------------------------------------------------- */
point **make_path(char *chr, point *s, point *t) {
    int i, x, y, z, disp, sign_x, sign_y, sign_z;

    // Define relative start and end points of coordinate values
    unsigned int x_begin = 0,           x_end = COORD_SIZE-1;
    unsigned int y_begin = x_end+1,     y_end = (2*(COORD_SIZE))-1;
    unsigned int z_begin = y_end+1,     z_end = (3*(COORD_SIZE))-1;

    point **path = (point**) malloc(sizeof(point*) * (N_POINTS+2));
    for(i = 0, disp = 0; i < N_POINTS+2; i++, disp+=(z_end+1))
        path[i] = (point*) malloc(sizeof(point));

    // Convert the string into a series of (x,y,z) coodinates
    pt_copy(path[0], s);
    for(i = 1, disp = 0; i < N_POINTS+1; i++, disp+=(z_end+1)) {
        sign_x = binToDecimal(chr, disp+x_begin, disp+x_begin);
        sign_y = binToDecimal(chr, disp+y_begin, disp+y_begin);
        sign_z = binToDecimal(chr, disp+z_begin, disp+z_begin);

        x = binToDecimal(chr, disp+x_begin+1, disp+x_end);
        y = binToDecimal(chr, disp+y_begin+1, disp+y_end);
        z = binToDecimal(chr, disp+z_begin+1, disp+z_end);

        if (sign_x) path[i]->x = x*-1;  else path[i]->x = x;
        if (sign_y) path[i]->y = y*-1;  else path[i]->y = y;
        if (sign_z) path[i]->z = z*-1;  else path[i]->z = z;
    }
    pt_copy(path[N_POINTS+1], t);

    return path;
}


void free_path(point **path) {
    int i;
    for (i = 0; i < N_POINTS+2; i++) free(path[i]);
    free(path);
}


/* ------------------------------------------------------------------------- */
/* Convert a substring of binary digits to a decimal value                   */
/* Uses the following recursive formula:                                     */
/* f(0) = 0                                                                  */
/* f(1) = 1                                                                  */
/* f(n0) = 2*bin(n)                                                          */
/* f(n1) = 2*bin(n) + 1                                                      */
/*                                                                           */
/* bin      : An array of characters representing a bit string               */
/* begin    : The begin index (inclusive) of the substring to convert        */
/* end      : The end index (inclusive) of the substring to convert          */
/* ------------------------------------------------------------------------- */
int binToDecimal(char* bin, unsigned int begin, unsigned int end) {
    if (begin > end || begin < 0) return 0;
    int i, digit, dec_val = 0;

    for (i = begin; i <= end; i++) {
        digit = (int) (bin[i] - '0');
        dec_val = 2*dec_val + digit;
    }

    return dec_val;
}

