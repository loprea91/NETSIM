#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

/* defined constants */
#define FILE_BUFFER_SIZE 2048

enum NETSIM_ERROR_CODES {
    NETSIM_NOERROR = 0,
    NETSIM_ERROR = 1,
    NETSIM_BAD_PARAMETER = 2,
    NETSIM_EARLY_EXIT = -1
};

/* individual parameter structures */
struct neuron_parameters {
    double taum;
    double taur;
    double tauref;
    double vr;
    double vth;
    double vreset;
    double Gl;
    double El;
    double Ie;
    double Cm;
};

struct synapse_parameters {
    double ge;
    double gi;
    double taue;
    double taui;
    double Ee;
    double Ei;
    double synapse_delay;
    double min_conduction_speed;
    double max_conduction_speed;
    double p_release;
    double min_uniform_delay;
    double max_uniform_delay;
};

struct network_parameters {
    unsigned int N;
    unsigned int K;
    unsigned int Ne;
    unsigned int Ni;
    unsigned int Nepop;
    unsigned int Kin;
    unsigned int Kout;
    unsigned int spatial_dimensions;
    unsigned int poisson_K_E;
    unsigned int poisson_K_I;

    double L;
    double sigma_space;
    double poisson_rate;
    double poisson_start_time;
    double poisson_stop_time;
    double rewiring_probability;

    unsigned int NrowE;
    unsigned int NrowI;
    double dx_exc;
    double dx_inh;
};

struct simulation_parameters {
    /* simulation control parameters */
    int connector;
    int initiator;
    int conduction_speed_code;
   
    double T;
    double dt;
	
    double vm_mean;
    double vm_sigma;
    double ge_mean;
    double ge_sigma;
    double gi_mean;
    double gi_sigma;
	
    int gamma_parameter_shape;
    double gamma_parameter_scale;
	
    char *parameter_file_path;
    char *optional_connection_path;
    char *output_path;
    char *current_injection_path;
    char *sparse_current_injection_path;
    char *poisson_rate_path;

    int buffer_length;

    int job_id;
	
    /* external inputs to simulation */
    int current_injection_count;
    double *current_injection;
    uint32_t sparse_current_entry_count;
    struct sparse_current_injection_entry *sparse_injection;
    double *poisson_variable_rate;

    /* single neuron stimulation parameters (DS3) */
    bool stim;
    int stim_id;
    double stim_time;
    double perturbation;


    /* recording parameters */
    int output_path_override;
    int connection_seed;
    int output_file_code;
    int save_connectivity;
    int save_positions;
    int record_downsample_factor;
    double start_record_time;
    double stop_record_time;
    int partial_recording_time;
};

/* umbrella parameter structure */
struct parameters {
    struct neuron_parameters neuron;
    struct synapse_parameters synapse;
    struct network_parameters network;
    struct simulation_parameters simulation;
};
    
/* neuron structure */
struct neuron {
    bool *excitatory;
    double *v;
    double *ge;
    double *gi;
    uint32_t *K;
    double *conduction_speed;
    double *x_position;
    double *y_position;
    double *z_position;
    double *lastSpike;
    int *lastPoisE;
    int *lastPoisI;
    double *buff_input_E;
    double *buff_input_I;
};

/* mappings for spare current injection format*/
struct sparse_current_injection_entry {
    uint32_t start_step;
    uint32_t stop_step;
    uint32_t index_count;
    uint32_t *indices;
    uint32_t pos;
    double intensity;
};

#endif
