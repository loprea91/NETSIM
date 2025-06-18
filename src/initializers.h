#ifndef _INITIALIZERS_H_
#define _INITIALIZERS_H_

#include "simulation.h"

enum NETSIM_INIT_TYPE { NETSIM_INIT_SUSTAINED = 0, NETSIM_INIT_NONE = 1, NETSIM_INIT_POISSON = 2 };

enum NETSIM_CONDITION_SPEED_CODE {
    NETSIM_CONSTANT_CODE = 0,
    NETSIM_UNIFORM_CODE = 1,
    NETSIM_GAMMA_CODE = 2,
};

int initialize_network(struct parameters *p, struct neuron *network, struct rng_state *rng);
int allocate_input_buffers(struct parameters *p, struct neuron *network);

int initialize_sustained_simulation(struct parameters *p, struct neuron *network,
                                    struct rng_state *rng);

int initialize_1D_ring_simulation(struct parameters *p, struct neuron *network);
int initialize_2D_lattice_simulation(struct parameters *p, struct neuron *network);

#endif
