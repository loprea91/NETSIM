#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rng.h"
#include "rng_dist.h"

#include "connectors.h"
#include "initializers.h"

int connect_network(struct parameters *p, struct neuron *network, struct rng_state *rng) {
    /* make connections */
    printf("     Connecting network...\n");
    fflush(stdout);

    if (p->simulation.connector == NETSIM_RANDOM_CONNECT) {
        initialize_1D_ring_simulation(p, network);
        random_connect(&p->network, &p->simulation, &p->synapse, network, rng);
        p->network.spatial_dimensions = 1;
    } else if (p->simulation.connector == NETSIM_GAUSSIAN_CONNECT) {
        initialize_1D_ring_simulation(p, network);
        gaussian_connect(&p->network, &p->simulation, &p->synapse, network, rng);
        p->network.spatial_dimensions = 1;
    } else if (p->simulation.connector == NETSIM_TEST_CONNECT) {
        test_connect(p, network);
    } else if (p->simulation.connector == NETSIM_GAUSSIAN_2D_CONNECT) {
        initialize_2D_lattice_simulation(p, network);
        gaussian_connect_2D(&p->network, &p->simulation, &p->synapse, network, rng);
        p->network.spatial_dimensions = 2;
    } else if (p->simulation.connector == NETSIM_CLUSTERED_CONNECT) {
        initialize_1D_ring_simulation(p, network);
        // initialize_2D_lattice_simulation( p, network );
        clustered_connect(&p->network, &p->simulation, &p->synapse, network, rng);
        p->network.spatial_dimensions = 2;
    } else {
        printf("Given connector not specified\n");
        return NETSIM_BAD_PARAMETER;
    }

    /* random rewiring if specified */
    if (p->network.rewiring_probability > 0) {
        printf("     Rewiring network with probability %2.2f...\n",
               p->network.rewiring_probability);
        fflush(stdout);
        random_rewiring(&p->network, &p->simulation, &p->synapse, network, rng);
    }

    return NETSIM_NOERROR;
}

// Calculates shortest distance bewteen two points on a ring with length L
// if distance bewteen each point is equally spaced, dx
double calculate_shortest(int i, int j, double dxs, double dxr, double L) {
    double opt1 = fabs(i * dxs - j * dxr);
    double opt2 = (L - opt1);
    return ((opt1 < opt2) ? opt1 : opt2);
}

/* function: calculate shortest distance between two nodes on a circle. */
/* warning: distance given has satisfy  distance/dxj > -N */
int calculate_closest(double distance, double position, double dxj, int N) {
    double jdis;
    int j;

    jdis = position + distance;
    j = round(jdis / dxj);

    j = (j + N) % N;

    return j;
}

/* gaussian connect */
int gaussian_connect(struct network_parameters *net, struct simulation_parameters *sim,
                     struct synapse_parameters *syn, struct neuron *network,
                     struct rng_state *rng) {
    /* init */
    const int N = net->N;
    const int Ne = net->Ne;
    const int Ni = net->Ni;
    const int K = net->K;
    const double dx_exc = net->dx_exc;
    const double dx_inh = net->dx_inh;

    int count, count_connections;
    double x_position, distance;
    double max_distance = 0;

    /* temporarily save connections in boolean size N, then translate to sparse */
    /* like this it is way faster to check for double connections */
    bool temparr[N];

    int jj = 0;
    for (int ii = 0; ii < net->N; ii++) {
        x_position = network->x_position[ii];

        /* reset temparr */
        memset(temparr, (bool)0, (size_t)N * sizeof(bool));

        count = (int)round(.8 * K);
        count_connections = 0;
        while (count > 0) {
            /* draw gaussian distance */
            distance = rng_gauss(rng) * net->sigma_space;

            /* within range? */
            if (fabs(distance) > (.5 * net->L)) {
                continue;
            }

            jj = calculate_closest(distance, x_position, dx_exc, Ne);

            /* no self- or double connections */
            if (((ii == jj) && network->excitatory[ii]) || temparr[jj]) {
                continue;
            }

            /* check for the largest connection distance */
            if (fabs(distance) > max_distance) {
                max_distance = fabs(distance);
            }

            /* connect */
            network->K[(ii * K) + count_connections] = jj;

            temparr[jj] = true;
            count--;
            count_connections++;
        }

        count = (int)round(.2 * net->K);
        while (count > 0) {
            distance = rng_gauss(rng) * net->sigma_space;

            if (fabs(distance) > (.5 * net->L)) {
                continue;
            }

            jj = Ne + calculate_closest(distance, x_position, dx_inh, Ni);

            /* no self- or double connections */
            if (((ii == jj) && !network->excitatory[ii]) || temparr[jj]) {
                continue;
            }

            /* check for the largest connection distance */
            if (fabs(distance) > max_distance) {
                max_distance = fabs(distance);
            }

            /* connect */
            network->K[(ii * K) + count_connections] = jj;

            temparr[jj] = true;
            count--;
            count_connections++;
        }
    }

    /* assign buffer length based on max distance */
    sim->buffer_length =
        (int)round((syn->synapse_delay + (max_distance / (syn->min_conduction_speed))) / sim->dt) +
        1;

    return NETSIM_NOERROR;
}

/* gaussian connect */
int gaussian_connect_2D(struct network_parameters *net, struct simulation_parameters *sim,
                        struct synapse_parameters *syn, struct neuron *network,
                        struct rng_state *rng) {
    /* init */
    int ii, jj, c1, c2;
    int count, count_connections;
    double x_position, y_position, distance_x, distance_y, distance;

    double max_distance = 0;

    /* calculate */
    const int Ne = net->Ne;
    const int K = net->K;
    const int NrowE = net->NrowE;
    const int NrowI = net->NrowI;
    const double dx_exc = net->dx_exc;
    const double dx_inh = net->dx_inh;

    /* temporarily save connections in boolean size N, then translate to sparse */
    /* like this it is way faster to check for double connections */
    bool temparr[net->N];

    for (ii = 0; ii < net->N; ii++) {
        x_position = network->x_position[ii];
        y_position = network->y_position[ii];

        /* reset temparr */
        memset(temparr, (bool)0, (size_t)net->N * sizeof(bool));

        count = (int)round(.8 * net->K);
        count_connections = 0;
        while (count > 0) {

            /* draw gaussian distance */
            distance_x = rng_gauss(rng) * net->sigma_space;
            distance_y = rng_gauss(rng) * net->sigma_space;
            distance = sqrt(pow(distance_x, 2) + pow(distance_y, 2));

            /* within range? */
            if (fabs(distance) > (M_SQRT2 * .5 * net->L)) {
                continue;
            }

            c1 = calculate_closest(distance_x, x_position, dx_exc, NrowE); /* x-coordinate */
            c2 = calculate_closest(distance_y, y_position, dx_exc, NrowE); /* y-coordinate */
            jj = c1 + (c2 * NrowE);

            /* no self- or double connections */
            if (((ii == jj) && network->excitatory[ii]) || temparr[jj]) {
                continue;
            }

            /* check for the largest connection distance */
            if (fabs(distance) > max_distance) {
                max_distance = fabs(distance);
            }

            /* connect */
            network->K[(ii * K) + count_connections] = jj;

            temparr[jj] = true;
            count--;
            count_connections++;
        }

        count = (int)round(.2 * net->K);
        while (count > 0) {
            /* draw gaussian distance */
            distance_x = rng_gauss(rng) * net->sigma_space;
            distance_y = rng_gauss(rng) * net->sigma_space;
            distance = sqrt(pow(distance_x, 2) + pow(distance_y, 2));

            if (fabs(distance) > (M_SQRT2 * .5 * net->L)) {
                continue;
            }

            c1 = calculate_closest(distance_x, x_position, dx_inh, NrowI); /* x-coordinate */
            c2 = calculate_closest(distance_y, y_position, dx_inh, NrowI); /* y-coordinate */
            jj = Ne + c1 + (c2 * NrowI);

            /* no self- or double connections */
            if (((ii == jj) && !network->excitatory[ii]) || temparr[jj]) {
                continue;
            }

            /* check for the largest connection distance */
            if (fabs(distance) > max_distance) {
                max_distance = fabs(distance);
            }

            /* connect */
            network->K[(ii * K) + count_connections] = jj;

            temparr[jj] = true;
            count--;
            count_connections++;
        }
    }

    /* assign buffer length based on max distance */
    sim->buffer_length =
        (int)round((syn->synapse_delay + (max_distance / (syn->min_conduction_speed))) / sim->dt) +
        1;

    /* cleanup */
    return NETSIM_NOERROR;
}

/* function: Erdos-Renyi random graph */
int random_connect(struct network_parameters *net, struct simulation_parameters *sim,
                   struct synapse_parameters *syn, struct neuron *network, struct rng_state *rng) {
    /* init */
    const int N = net->N;
    const int K = net->K;
    uint32_t target = 0;

    /* loop over network */
    for (int ii = 0; ii < N; ii++) {
        for (int jj = 0; jj < K; jj++) {
            /* connect */
            target = rng_uint32(rng) % N;
            network->K[(ii * K) + jj] = target;
        }
    }

    sim->buffer_length =
        (int)((syn->synapse_delay + ((.5 * net->L) / (syn->min_conduction_speed))) / sim->dt) + 1;

    return NETSIM_NOERROR;
}

/* function: clustered connectivity graph */
int clustered_connect(struct network_parameters *net, struct simulation_parameters *sim,
                      struct synapse_parameters *syn, struct neuron *network,
                      struct rng_state *rng) {
    /* init */
    const int N = net->N;
    const int K = net->K;
    const int Kin = net->Kin;
    const int Kout = net->Kout;
    const int Nepop = net->Nepop;
    const int Ne = net->Ne;

    /* loop over excitatory population */
    for (int ii = 0; ii < Ne; ii++) {
        const int cluster_ii = floor(ii / Nepop);

        for (int jj = 0; jj < Kin; jj++) {
            const int target = (rng_uint32(rng) % Nepop) + (cluster_ii * Nepop);
            network->K[(ii * K) + jj] = target;
        }

        for (int jj = Kin; jj < Kin + Kout; jj++) {
            const int target = (rng_uint32(rng) % net->N);
            network->K[(ii * K) + jj] = target;
        }

        for (int jj = Kin + Kout; jj < K; jj++) {
            const int target = (rng_uint32(rng) % (N - Ne)) + Ne;
            network->K[(ii * K) + jj] = target;
        }
    }

    /* loop over inhibitory population */
    for (int ii = Ne; ii < net->N; ii++) {
        for (int jj = 0; jj < net->K; jj++) {
            const int target = rng_uint32(rng) % net->N;
            network->K[(ii * K) + jj] = target;
        }
    }

    return NETSIM_NOERROR;
}

/* function: load test conections from file */
int test_connect(struct parameters *p, struct neuron *network) {
    /* init */
    int ii, jj, id;
    int N, K;
    FILE *file;

    N = p->network.N;
    K = p->network.K;

    file = fopen(p->simulation.optional_connection_path, "r");
    int r = 0;

    /* loop over network */
    for (ii = 0; ii < N; ii++) {
        for (jj = 0; jj < K; jj++) {
            /* get connection from file */
            r = fscanf(file, "%d\n", &id);
            if (r < 0) {
                // throw and error?
                // fclose( file );
                // return NETSIM_ERROR;
            }
            /* connect */
            network->K[ii * K + jj] = id;
        }
    }

    /* clean up */
    fclose(file);
    return NETSIM_NOERROR;
}

/* function: random rewiring of EXCITATORY network connections */
int random_rewiring(struct network_parameters *net, struct simulation_parameters *sim,
                    struct synapse_parameters *syn, struct neuron *network, struct rng_state *rng) {

    /* init */
    const int N = net->N;
    const int Ne = net->Ne;
    const int K = net->K;
    uint32_t target = 0;

    /* loop over network */
    for (int ii = 0; ii < Ne; ii++) {
        for (int jj = 0; jj < K; jj++) {
            if (rng_dbl64(rng) < net->rewiring_probability) {
                /* select random target */
                target = rng_uint32(rng) % N;

                /* write random target */
                network->K[(ii * K) + jj] = target;
            }
        }
    }

    /* clean up */
    return NETSIM_NOERROR;
}
