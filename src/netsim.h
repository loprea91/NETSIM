#ifndef _NETSIM_H_
#define _NETSIM_H_

#include "output.h"
#include "simulation.h"

int *restrict const get_distance_matrix(const int N, const int K, const int spatial_dimension,
					const double synapse_delay, const double L, const double dt,
					const uint32_t *restrict const K_conn,
					const double *restrict const conduction_speed,
					const double *restrict const x_position,
					const double *restrict const y_position);
int run_simulation(struct parameters *p, struct neuron *network, struct output *o,
                   struct rng_state *rng);

#endif
