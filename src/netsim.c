/********************************************************************************
 *                                                                              *
 *  FILE:   netsim.c                                                            *
 *  VERSION:  0.1                                                               *
 *  PROGRAM:  NETSIM                                                            *
 *                                                                              *
 *  PURPOSE:  A fast, large-scale simulator for topographic spiking networks    *
 *                                                                              *
 *  Copyright (C) 2016-2020 Lyle Muller                                         *
 *  http://mullerlab.ca                                                         *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 *  DEVELOPMENT: Lyle Muller, Charlee Fletterman, Theo Desbordes, Gabriel       *
 *  Benigno, Christopher Steward                                                *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 * This file is part of NETSIM.                                                 *
 *                                                                              *
 *     NETSIM is free software: you can redistribute it and/or modify           *
 *     it under the terms of the GNU General Public License as published by     *
 *     the Free Software Foundation, either version 3 of the License, or        *
 *     (at your option) any later version.                                      *
 *                                                                              *
 *     NETSIM is distributed in the hope that it will be useful,                *
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *     GNU General Public License for more details.                             *
 *                                                                              *
 *     You should have received a copy of the GNU General Public License        *
 *     along with NETSIM.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                              *
 ********************************************************************************/

/* global includes */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/* local includes */
#include "isaac64.h"
#include "rng.h"
#include "rng_dist.h"

#include "netsim.h"

int *restrict const get_distance_matrix(const int N, const int K, const int spatial_dimensions,
					const double synapse_delay, const double L, const double dt,
					const uint32_t *restrict const K_conn,
					const double *restrict const conduction_speed,
					const double *restrict const x_position,
					const double *restrict const y_position)
{
    int *restrict const distance_matrix = calloc(N * K, sizeof(int));
    if(spatial_dimensions == 1) {
	for (int ii = 0; ii < N; ii++) {
	    for (int jj = 0; jj < K; jj++) {
		const int target = K_conn[(ii * K) + jj];
		const double diff_x = fabs(x_position[ii] - x_position[target]);
		const double distance_x = fmin(L - diff_x, diff_x);
		const double distance = sqrt((distance_x * distance_x));
		const int delay_in_bins =
		    (int)round((synapse_delay + (distance / conduction_speed[ii])) / dt);
		distance_matrix[(ii * K) + jj] = delay_in_bins;
	    }
	}
    }
    else if(spatial_dimensions == 2) {
	for (int ii = 0; ii < N; ii++) {
	    for (int jj = 0; jj < K; jj++) {
		const int target = K_conn[(ii * K) + jj];
		const double diff_x = fabs(x_position[ii] - x_position[target]);
		const double distance_x = fmin(L - diff_x, diff_x);
		const double diff_y = fabs(y_position[ii] - y_position[target]);
		const double distance_y = fmin(L - diff_y, diff_y);
		const double distance = sqrt((distance_x * distance_x) + (distance_y * distance_y));
		const int delay_in_bins =
		    (int)round((synapse_delay + (distance / conduction_speed[ii])) / dt);
		distance_matrix[(ii * K) + jj] = delay_in_bins;
	    }
	}
    }
    else {
	printf("Spacial dimension set incorrectly!\n");
    }
    
    return distance_matrix;
}

/* function: run_simulation */
int run_simulation(struct parameters *p, struct neuron *network, struct output *o,
                   struct rng_state *rng) {
    /* make local structure copies */
    const struct network_parameters net = p->network;
    const struct neuron_parameters neu = p->neuron;
    const struct synapse_parameters syn = p->synapse;
    const struct simulation_parameters sim = p->simulation;

    double *restrict const v = network->v;
    double *restrict const ge = network->ge;
    double *restrict const gi = network->gi;
    double *restrict const buff_input_E = network->buff_input_E;
    double *restrict const buff_input_I = network->buff_input_I;
    const uint32_t *restrict const K_conn = network->K;
    const bool *restrict const excitatory = network->excitatory;
    const double *restrict const x_position = network->x_position;
    const double *restrict const y_position = network->y_position;
    const double *restrict const conduction_speed = network->conduction_speed;
    double *restrict const lastSpike = network->lastSpike;

    /* init local variables */
    const double T = sim.T;
    const double dt = sim.dt;
    const double ndt = -sim.dt;
    const int buffer_length = sim.buffer_length - 1;

    const unsigned int N = net.N;
    const unsigned int Nsteps = round(T/dt);
    const unsigned int K = net.K;
    const double L = net.L;

    const double we = syn.ge;
    const double wi = syn.gi;
    const double Ee = syn.Ee;
    const double Ei = syn.Ei;
    const double taue = syn.taue;
    const double taui = syn.taui;
    const double synapse_delay = syn.synapse_delay;

    const double taur = neu.taur;
    const double taum = neu.taum;
    const double El = neu.El;
    const double Cm = neu.Cm;
    const double Gl = Cm/taum;
    const double vreset = neu.vreset;
    const double vth = neu.vth;

#if defined(CURRENT_INJECTION) || defined(SPARSE_CURRENT_INJECTION)
    double Ie = neu.Ie;
#else
    const double Ie = 0;
#endif

#if defined(EXTERNAL_INPUT)
    const double poisson_rate = net.poisson_rate;
#endif
#if defined(VARIABLE_EXTERNAL_INPUT)
    double poisson_rate = net.poisson_rate;
     const double *poisson_variable_rate = sim.poisson_variable_rate;
#endif
    
#if defined(EXTERNAL_INPUT) || defined(VARIABLE_EXTERNAL_INPUT)
    const unsigned int poisson_start_time_steps = (unsigned int)round(net.poisson_start_time/dt);
    const unsigned int poisson_stop_time_steps = (unsigned int)round(net.poisson_stop_time/dt);
    const unsigned int poisson_K_E = net.poisson_K_E;
    const unsigned int poisson_K_I = net.poisson_K_I;
#endif

#if defined(RELEASE_PROBABILITY)
    const double p_release = syn.p_release;
    unsigned int spike_transmitted = 0;
    unsigned int spike_not_transmitted = 0;
#endif

#if defined(CURRENT_INJECTION)
    const double *current_injection = sim.current_injection;
    const int current_injection_map = N/sim.current_injection_count;
#endif
    
    int *restrict const distance_matrix =
	get_distance_matrix(N, K, p->network.spatial_dimensions,
			    synapse_delay, L, dt, 
			    K_conn, conduction_speed,
			    x_position, y_position);
    
    /* init recording */
    struct recorder rec;
#if defined(LONG_SIMULATION)
    initialize_partial_recorders(&rec, p);
#else
    initialize_recorders(&rec, p);
#endif
    
    const int record_downsample_factor = rec.downsample_factor;
    const int record_start_bin = rec.start_bin;
    const int record_stop_bin = rec.stop_bin;
    int print_steps = (int)(1.0/dt);
    /* MAIN SIMULATION LOOP */
    for (int tt = 0; tt < Nsteps; tt++) {
        const double t = dt * tt;
        const int buffer_position = tt & buffer_length;
        const int buffer_idx = buffer_position * N;

        /* Print updates */
        double elapsed = tt * dt;
        double progress = ((double)tt / Nsteps) * 100.0;
        if (tt % print_steps == 0) {
            printf("Simulation time: %.2f seconds, Progress: %.2f%%\n", elapsed, progress);
            fflush(stdout);
        }

        /* add external inputs to ge or gi buffers */
        for (int ii = 0; ii < N; ii++) {
#include "templates/external_input.h"
#include "templates/variable_external_input.h"
        }

        /* integrate synaptic input and then reset buffer */
        for (int ii = 0; ii < N; ii++) {
	    ge[ii] += ndt * (ge[ii] / taue) + buff_input_E[buffer_idx + ii];
	    buff_input_E[buffer_idx + ii] = 0;
	    
            gi[ii] += ndt * (gi[ii] / taui) + buff_input_I[buffer_idx + ii];
	    buff_input_I[buffer_idx + ii] = 0;
	}
	
        for (int ii = 0; ii < N; ii++) {
            if (t > (lastSpike[ii] + taur)) {
/* setup Ie if there is current injection*/
#include "templates/current_injection.h"
#include "templates/sparse_current_injection.h"

                /* integrate membrane equation */
                v[ii] += dt * ((ge[ii] * (Ee - v[ii]) + /* excitatory synapse */
                                gi[ii] * (Ei - v[ii]) + /* inhibitory synapse */
                                Gl     * (El - v[ii]) + /* leak conductance */
                                Ie) /
                               Cm);
            }
        }

        for (int ii = 0; ii < N; ii++) {
            if (v[ii] >= vth) {
		    
                /* communicate synaptic outputs */
                for (int jj = 0; jj < K; jj++) {
#include "templates/release_probability_template_start.h"
                    const int target = K_conn[(ii * K) + jj];
                    const int delay_in_bins = distance_matrix[(ii * K) + jj];
                    const int insert_position = (buffer_position + delay_in_bins) & buffer_length;
                    const int buff_input_idx = (insert_position * N) + target;
#include "templates/release_probability_template_stop.h"

                    if (excitatory[ii]) {
                        buff_input_E[buff_input_idx] += we;
                    }
		    else {
                        buff_input_I[buff_input_idx] += wi;
                    }
                }
		
                /* record spikes */
                record_spike(&rec, N, T, t, ii);
		
                /* reset vm and record spike time */
                v[ii] = vreset;
                lastSpike[ii] = t;
            }
        }
	
	/* recording */
	record_state_variables(&sim, &rec, network, N, tt);
    }

    /* record spike count to output */
    o->spike_count = rec.spike_count;

    /* collect static variables */
    record_static_variables(&rec, network, N);

    /* write recorded values to file */
    write_output_files(&sim, &rec, N);
    
    /* free distance matrix */
    free(distance_matrix);

    destory_recorders(&rec, p);

    /* return value */
    return NETSIM_NOERROR;
}
