#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "isaac64.h"
#include "rng.h"
#include "rng_dist.h"

#include "initializers.h"

int initialize_network(struct parameters *p, struct neuron *network, struct rng_state *rng) {
    /* init */
    int ii;
    int N = p->network.N;
    int Ne = (int)round(.8 * N);
    int Ni = p->network.N - Ne;
    double tmp_conduction_speed;

    p->network.Ne = Ne;
    p->network.Ni = Ni;

    /* allocate space */
    network->excitatory = calloc(N, sizeof(bool));
    network->v = calloc(N, sizeof(double));
    network->ge = calloc(N, sizeof(double));
    network->gi = calloc(N, sizeof(double));
    network->conduction_speed = calloc(N, sizeof(double));
    network->x_position = calloc(N, sizeof(double));
    network->y_position = calloc(N, sizeof(double));
    network->z_position = calloc(N, sizeof(double));
    network->lastSpike = calloc(N, sizeof(double));

    assert(network->excitatory != NULL &&
	   network->v != NULL &&
	   network->ge != NULL &&
           network->gi != NULL &&
	   network->conduction_speed != NULL &&
           network->x_position != NULL &&
	   network->y_position != NULL &&
           network->z_position != NULL &&
	   network->lastSpike != NULL);

    /* initialize network */
    for (ii = 0; ii < p->network.N; ii++) {
        network->excitatory[ii] = (ii < Ne);
        network->v[ii] = p->neuron.El;
        network->ge[ii] = 0;
        network->gi[ii] = 0;
        network->lastSpike[ii] = -1;

        /* calculate conduction velocity: code 0 -> constant speed, 1 -> uniform, 2 -> gamma
         * distribution  3 -> multilayer*/
        if (p->simulation.conduction_speed_code == NETSIM_CONSTANT_CODE) {
            network->conduction_speed[ii] = p->synapse.min_conduction_speed;
        } else if (p->simulation.conduction_speed_code == NETSIM_UNIFORM_CODE) {
            network->conduction_speed[ii] =
                (p->synapse.max_conduction_speed - p->synapse.min_conduction_speed) *
                    rng_dbl64(rng) +
                p->synapse.min_conduction_speed;
        } else if (p->simulation.conduction_speed_code == NETSIM_GAMMA_CODE) {
            do {
                tmp_conduction_speed = gamma_sampling(p->simulation.gamma_parameter_shape,
                                                      p->simulation.gamma_parameter_scale, rng);
            } while ((tmp_conduction_speed < p->synapse.min_conduction_speed) ||
                     (tmp_conduction_speed > p->synapse.max_conduction_speed));

            network->conduction_speed[ii] = tmp_conduction_speed;
        }
    }

    return NETSIM_NOERROR;
}

/* function: allocate input buffers */
int allocate_input_buffers(struct parameters *p, struct neuron *network) {

    /* init */
    int buffer_length;

    /* calculate needed buffer size */
    if (p->synapse.max_uniform_delay > 0) {
        buffer_length = (int)(p->synapse.max_uniform_delay / p->simulation.dt) + 1;
        p->simulation.buffer_length = buffer_length;
    }

    p->simulation.buffer_length = 1 << (64 - __builtin_clzl(p->simulation.buffer_length - 1));

    buffer_length = p->simulation.buffer_length;

    /* initialize buffers in network */
    network->buff_input_E = calloc(p->network.N * buffer_length, sizeof(double));
    network->buff_input_I = calloc(p->network.N * buffer_length, sizeof(double));

    /* check memory allocation */
    assert(network->buff_input_E != NULL &&
	   network->buff_input_I != NULL);

    /* initialize network */
    for (int ii = 0; ii < buffer_length; ii++) {
        network->buff_input_E[ii] = 0;
        network->buff_input_I[ii] = 0;
    }

    return NETSIM_NOERROR;
}

/* function: intialize for self-sustained activity */
int initialize_sustained_simulation(struct parameters *p, struct neuron *network,
                                    struct rng_state *rng) {
    printf("> Initializing self-sustained activity ");

    if (p->simulation.initiator == NETSIM_INIT_SUSTAINED) {
        printf("with gaussian vm, ge, and gi\n");
        /* init */
        for (int ii = 0; ii < p->network.N; ii++) {
            /* init v, ge and gi with normal distribution for sustained (emperical) */
            network->v[ii] = rng_gauss(rng) * p->simulation.vm_sigma + p->simulation.vm_mean;
            network->ge[ii] = rng_gauss(rng) * p->simulation.ge_sigma + p->simulation.ge_mean;
            network->gi[ii] = rng_gauss(rng) * p->simulation.gi_sigma + p->simulation.gi_mean;
        }
    } else if (p->simulation.initiator == NETSIM_INIT_POISSON) {
    }

    return 0;
}

/* function: init to zero */
int initialize_1D_ring_simulation(struct parameters *p, struct neuron *network) {
    /* init */
    int N = p->network.N;
    int Ne = p->network.Ne;
    int Ni = p->network.Ni;
    double dx_exc, dx_inh;

    network->K = calloc(p->network.N * p->network.K, sizeof(uint32_t));
    assert(network->K != NULL);

    /* calculate */
    dx_exc = p->network.L / Ne;
    dx_inh = p->network.L / Ni;

    p->network.dx_exc = dx_exc;
    p->network.dx_inh = dx_inh;

    for (int ii = 0; ii < N; ii++) {
        /* calculate neuron's x_position in the ring */
        if (network->excitatory[ii] == true) {
            network->x_position[ii] = ii * dx_exc;
        } else {
            network->x_position[ii] = (ii - Ne) * dx_inh;
        }

        network->y_position[ii] = 0;
        network->z_position[ii] = 0;
    }

    return NETSIM_NOERROR;
}

int initialize_2D_lattice_simulation(struct parameters *p, struct neuron *network)
{
    /* init */
    const int Ne = p->network.Ne;
    const int Ni = p->network.Ni;

    network->K = calloc(p->network.N * p->network.K, sizeof(uint32_t));
    assert(network->K != NULL);

    /* calculate */
    const int NrowE = floor(sqrt(Ne));
    const int NrowI = floor(sqrt(Ni));

    const double dx_exc = p->network.L / NrowE;
    const double dx_inh = p->network.L / NrowI;

    p->network.NrowE = NrowE;
    p->network.NrowI = NrowI;

    p->network.dx_exc = dx_exc;
    p->network.dx_inh = dx_inh;

    /* loop over neurons and calculate lattice indices */
    for (int ii = 0; ii < NrowE; ii++) {
        for (int jj = 0; jj < NrowE; jj++) {
            /* calculate index */
            const int index = jj + (ii * NrowE);

            network->x_position[index] = jj * dx_exc;
            network->y_position[index] = ii * dx_exc;
        }
    }

    /* loop over neurons and calculate lattice indices */
    for (int ii = 0; ii < NrowI; ii++) {
        for (int jj = 0; jj < NrowI; jj++) {
            /* calculate index */
            const int index = jj + (ii * NrowI) + Ne;

            network->x_position[index] = jj * dx_inh;
            network->y_position[index] = ii * dx_inh;
        }
    }

    return NETSIM_NOERROR;
}
