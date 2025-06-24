#ifndef _CONNECTORS_H_
#define _CONNECTORS_H_

#include "simulation.h"

enum NETSIM_CONNECT_TYPE {
    NETSIM_RANDOM_CONNECT = 1,
    NETSIM_GAUSSIAN_CONNECT = 2,
    NETSIM_TEST_CONNECT = 3,
    NETSIM_GAUSSIAN_2D_CONNECT = 4,
    NETSIM_GAUSSIAN_MULTILAYER_2D_CONNECT = 5,
    NETSIM_RANDOM_HORIZTONAL_MULTILAYER_CONNECT = 6,
    NETSIM_CLUSTERED_CONNECT = 7
};

int connect_network(struct parameters *p, struct neuron *network, struct rng_state *rng);

double calculate_shortest(int i, int j, double dxs, double dxr, double L);
int calculate_closest(double distance, double position, double dxj, int N);
int gaussian_connect(struct network_parameters *net, struct simulation_parameters *sim,
                     struct synapse_parameters *syn, struct neuron *network, struct rng_state *rng);
int gaussian_connect_2D(struct network_parameters *net, struct simulation_parameters *sim,
                        struct synapse_parameters *syn, struct neuron *network,
                        struct rng_state *rng);
int random_connect(struct network_parameters *net, struct simulation_parameters *sim,
                   struct synapse_parameters *syn, struct neuron *network, struct rng_state *rng);
int test_connect(struct parameters *p, struct neuron *network);
int random_rewiring(struct network_parameters *net, struct simulation_parameters *sim,
                    struct synapse_parameters *syn, struct neuron *network, struct rng_state *rng);
int clustered_connect(struct network_parameters *net, struct simulation_parameters *sim,
                      struct synapse_parameters *syn, struct neuron *network,
                      struct rng_state *rng);

#endif
