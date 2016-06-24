#include "advect.h"
#define DIM_MAX (3)

// update particle position using 4th order Runge-Kutta method
void update_particle_position_rk4( e *E, double* cart_pos, double* vel, double dt, int dim){

        int d;
        double accumPos[DIM_MAX];       // where the pos[] solution vector is accumulated
        double pos[DIM_MAX];  // position in longitude, latitude
        double posPred[DIM_MAX];  // position in meters
        double k1[DIM_MAX];
        int whereILive;
        double fac;

        fac = 1.0/6.0;

        for( d=0; d<dim; d++ ) {

                accumPos[d] = cart_pos[d]; // initial step

                k1[d] = dt * vel[d];
                posPred[d] = cart_pos[d] + 0.5 * k1[d];

                accumPos[d] = accumPos[d] + fac * k1[d]; // first step ( + 1/6 k1 )

        }

        // convert posPred back to lon lat so we can call get_owner_element
  			pos[1] = r2d(posPred[1]/R);	// calc lat first
  			pos[0] = r2d(posPred[0] / (R*cos(d2r(pos[1]))));	// now lon
        // find out which element containts this position
        whereILive = get_owner_element(E, pos);
        // update weights
  			calculate_interpolation_weights(&E->el[whereILive], E->xi, E->eta, pos);
  			// interpolate the nodal velocity to the current point
  			interpolate_point(&E->el[whereILive], vel);


        /* Reuse k1. Here k1 = k2 */
        for( d=0; d<dim; d++ ) {
                k1[d] = dt * vel[d];
                posPred[d] = cart_pos[d] + 0.5 * k1[d];
                accumPos[d] = accumPos[d] + 2.0*fac*k1[d];      /* second step ( + 1/3 k2 ) */
        }


        /* Put positions back in posPred */
        //for( d=0; d<dim; d++ ) {
        //        posPred[d] = cart_pos[d] + 0.5 * k1[d];
        //}

        // convert cart_pos back to lon lat
  			pos[1] = r2d(posPred[1]/R);	// calc lat first
  			pos[0] = r2d(posPred[0] / (R*cos(d2r(pos[1]))));	// now lon
        // find out which element containts this position
        //printf("\t\t2: posPred[0] = %f, posPred[1] = %f\n", posPred[0], posPred[1]);
        //printf("\t\t2: pos[0] = %f, pos[1] = %f\n", pos[0], pos[1]);
        whereILive = get_owner_element(E, pos);
  			calculate_interpolation_weights(&E->el[whereILive], E->xi, E->eta, pos);
  			// interpolate the nodal velocity to the current point
  			interpolate_point(&E->el[whereILive], vel);

        /* Reuse k1. Here k1 = k3 */
        for( d=0; d<dim; d++ ) {
                k1[d] = dt * vel[d];
                posPred[d] = cart_pos[d] + k1[d];
                accumPos[d] = accumPos[d] + 2.0*fac*k1[d];      /* third step ( + 1/3 k3 ) */
        }



        /* Put positions back in posPred */
        //for( d=0; d<dim; d++ ) {
        //        posPred[d] = cart_pos[d] + k1[d];
        //}


        // convert cart_pos back to lon lat
  			pos[1] = r2d(posPred[1]/R);	// calc lat first
  			pos[0] = r2d(posPred[0] / (R*cos(d2r(pos[1]))));	// now lon
        // find out which element containts this position
        whereILive = get_owner_element(E, pos);
  			calculate_interpolation_weights(&E->el[whereILive], E->xi, E->eta, pos);
  			// interpolate the nodal velocity to the current point
  			interpolate_point(&E->el[whereILive], vel);

        /* Reuse k1. Here k1 = k4 */
        for( d=0; d<dim; d++ ) {
                k1[d] = dt * vel[d];

                accumPos[d] = accumPos[d] + fac*k1[d];  /* fourth step ( + 1/6 k4 ) */

                // get new position
                cart_pos[d] = accumPos[d];

        }

        //printf("done\n"); fflush(stdout);
}
