#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "simulation.h"

struct recorder {
    bool record_flag;

    float *vm;
    float *ge;
    float *gi;

    double *position;
    double *spike_times;
    uint32_t *spike_ids;

    int buffer_pos;
    int downsample_factor;
    int start_bin;
    int stop_bin;
    int bin_size;
    size_t spike_count;

    double partial_time;
    int partial_bins;
    int partial_spike_count;
    int partial_file_count;
    int *partial_bin_indices;
    int partial_file_max;

    FILE *spike_times_file;
    FILE *spike_ids_file;
    FILE *ge_file;
    FILE *gi_file;
    FILE *vm_file;
    FILE *position_file;
};

/* output structure */
struct output {
    size_t spike_count;
};


int save_connectivity(struct parameters *p, struct neuron *network);

int save_out_param_file(struct parameters *p, int spike_count);

int initialize_recorders(struct recorder *rec, struct parameters *p);
int initialize_partial_recorders(struct recorder *rec, struct parameters *p);
int destory_recorders(struct recorder *rec, struct parameters *p);

int record_spike(struct recorder *rec, const int N, const double T, const double t, const uint32_t ii);
int record_state_variables(const struct simulation_parameters *sim,
			   struct recorder *rec, struct neuron *network,
			   const unsigned int N, const unsigned int tt);
int record_static_variables(struct recorder *rec, struct neuron *network, const unsigned int N);

int write_output_files(const struct simulation_parameters *sim,
		       struct recorder *rec, const unsigned int N);

int write_partial_variables(const struct simulation_parameters *sim, struct recorder *rec,
			    const unsigned int tt, const unsigned int N);


#endif
