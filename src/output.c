#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "output.h"

int save_out_param_file(struct parameters *p, int spike_count)
{
    char str[2048];
    sprintf(str, "%s/%08dparams.txt", p->simulation.output_path, p->simulation.output_file_code);

    /* Open file*/
    FILE *file = NULL;
    file = fopen(str, "w");

    /* Run info */
    fprintf(file, "#Spike count: %d\n", spike_count);

    /* Parameters */
    fprintf(file, "\n#Neuron Parameters\n");
    fprintf(file,
            "taum = %e\n"
            "taur = %e\n"
            "tauref = %e\n"
            "vr = %e\n"
            "vth = %e\n"
            "vreset = %e\n"
            "Gl = %e\n"
            "El = %e\n"
            "Ie = %e\n"
            "Cm = %e \n",
            p->neuron.taum, p->neuron.taur, p->neuron.tauref, p->neuron.vr, p->neuron.vth,
            p->neuron.vreset, p->neuron.Gl, p->neuron.El, p->neuron.Ie, p->neuron.Cm);

    fprintf(file, "\n#Synapse Parameters\n");
    fprintf(file,
            "ge = %e\n"
            "gi = %e\n"
            "taue = %e\n"
            "taui = %e\n"
            "Ee = %e\n"
            "Ei = %e\n"
            "synapse_delay = %e\n"
            "min_conducation_speed = %e\n"
            "max_conduction_speed = %e\n"
            "p_release = %e\n"
            "min_uniform_delay = %e\n"
            "max_uniform_delay = %e\n",
            p->synapse.ge, p->synapse.gi, p->synapse.taue, p->synapse.taui, p->synapse.Ee,
            p->synapse.Ei, p->synapse.synapse_delay, p->synapse.min_conduction_speed,
            p->synapse.max_conduction_speed, p->synapse.p_release, p->synapse.min_uniform_delay,
            p->synapse.max_uniform_delay);

    fprintf(file, "\n#Network Parameters\n");
    fprintf(file,
            "N = %d\n"
	    "Ne = %d\n"
	    "Ni = %d\n"
	    "NrowE = %d\n"
	    "NrowI = %d\n"
            "K = %d\n"
            "Nepop = %d\n"
            "Kin = %d\n"
            "Kout = %d\n"
            "spatial_dimensions = %d\n"
            "poisson_K_E = %d\n"
            "poisson_K_I = %d\n"
            "L = %lf\n"
            "sigma_space = %lf\n"
            "poisson_rate = %lf\n"
            "poisson_start_time = %lf\n"
            "poisson_stop_time = %lf\n"
            "rewiring_probability = %lf\n",
            p->network.N, p->network.Ne, p->network.Ni, p->network.NrowE, p->network.NrowI,
	    p->network.K, p->network.Nepop, p->network.Kin, p->network.Kout,
            p->network.spatial_dimensions, p->network.poisson_K_E,
            p->network.poisson_K_I, p->network.L, p->network.sigma_space, p->network.poisson_rate,
            p->network.poisson_start_time, p->network.poisson_stop_time,
            p->network.rewiring_probability);

    fprintf(file, "\n#Simulation Parameters\n");
    fprintf(file,
            "connector = %d\n"
            "initiator = %d\n"
            "output_file_code = %d\n"
            "connection_seed = %d\n"
            "record_downsample_factor = %d\n"
            "buffer_length = %d\n"
            "conduction_speed_code = %d\n"
            "gamma_parameter_shape = %d\n"
            "output_path_override = %d\n"
            "T = %e\n"
            "dt = %e\n"
            "vm_mean = %e\n"
            "vm_sigma = %e\n"
            "ge_mean = %e\n"
            "ge_sigma = %e\n"
            "gi_mean = %e\n"
            "gi_sigma = %e\n"
            "gamma_parameter_scale = %e\n"
            "parameter_file_path = %s\n"
            "optional_connection_path = %s\n"
            "output_path = %s\n"
            "current_injection_path = %s\n"
            "sparse_current_injection_path = %s\n"
            "poisson_rate_path = %s\n"
            "current_injection_count = %d\n"
            "spare_current_entry_count = %d\n"
            "save_connectivity = %d\n"
            "save_positions = %d\n"
            "start_record_time = %e\n"
            "stop_record_time = %e\n"
            "job_id = %d\n",
            p->simulation.connector, p->simulation.initiator, p->simulation.output_file_code, p->simulation.connection_seed,
            p->simulation.record_downsample_factor, p->simulation.buffer_length,
            p->simulation.conduction_speed_code, p->simulation.gamma_parameter_shape,
            p->simulation.output_path_override, p->simulation.T, p->simulation.dt,
            p->simulation.vm_mean, p->simulation.vm_sigma, p->simulation.ge_mean,
            p->simulation.ge_sigma, p->simulation.gi_mean, p->simulation.gi_sigma,
            p->simulation.gamma_parameter_scale, p->simulation.parameter_file_path,
            p->simulation.optional_connection_path, p->simulation.output_path,
            p->simulation.current_injection_path, p->simulation.sparse_current_injection_path,
            p->simulation.poisson_rate_path, p->simulation.current_injection_count,
            p->simulation.sparse_current_entry_count,
            p->simulation.save_connectivity,
            p->simulation.save_positions, p->simulation.start_record_time,
            p->simulation.stop_record_time, p->simulation.job_id);

    /* Compile time constants */
    fprintf(file, "\n#Compile time constants\n");
#ifdef CURRENT_INJECTION
    fprintf(file, "#CURRENT_INJEcTION=yes\n");
#endif

#ifdef SPARSE_CURRENT_INJECTION
    fprintf(file, "#SPARSE_CURRENT_INJECTION=yes\n");
#endif

#ifdef EXTERNAL_INPUT
    fprintf(file, "#EXTERNAL_INPUT=yes\n");
#endif

#ifdef VARIABLE_EXTERNAL_INPUT
    fprintf(file, "#VARIABLE_EXTERNAL_INPUT=yes\n");
#endif

#ifdef RELEASE_PROBABILITY
    fprintf(file, "#RELEASE_PROBABILITY=yes\n");
#endif

#if defined(LONG_SIMULATION)
    fprintf(file, "#LONG_SIMULATION=yes\n");
#endif

    fclose(file);

    return NETSIM_NOERROR;
}

int save_connectivity(struct parameters *p, struct neuron *network)
{

    printf("Saving connectivity...\n");

    /* init */
    char str[2048];

    /* bitfield variables */
    uint32_t send[p->network.K];
    uint32_t receive[p->network.K];

    sprintf(str, "%s/%08dii.bin", p->simulation.output_path, p->simulation.output_file_code);
    printf("Saving sending connections in: %s\n", str);
    FILE *s = fopen(str, "wb");

    sprintf(str, "%s/%08djj.bin", p->simulation.output_path, p->simulation.output_file_code);
    printf("Saving receiving connections in: %s\n", str);
    FILE *r = fopen(str, "wb");

    for (int ii = 0; ii < p->network.N; ii++) {
        // write array of repeating ii of same size as receiving neurons
        for (int jj = 0; jj < p->network.K; jj++) {
            send[jj] = ii;
            receive[jj] = network->K[(ii * p->network.K) + jj];
        }

        fwrite(send, sizeof(uint32_t), p->network.K, s);
        fwrite(receive, sizeof(uint32_t), p->network.K, r);
    }

    fclose(s);
    fclose(r);

    return NETSIM_NOERROR;
}
    
int initialize_recorders(struct recorder *rec, struct parameters *p)
{
    /* set recording flag */
    rec->record_flag = false;
    rec->spike_count = 0;
    rec->buffer_pos = 0;
        
    size_t size = p->network.N * 100 * p->simulation.T;
    char output_path[256];
    
    sprintf(output_path, "%s/%08dspk_times.bin", p->simulation.output_path, p->simulation.output_file_code);
    rec->spike_times_file = fopen(output_path, "wb");
    assert(rec->spike_times_file != NULL);
    
    sprintf(output_path, "%s/%08dspk_ids.bin", p->simulation.output_path, p->simulation.output_file_code);
    rec->spike_ids_file = fopen(output_path, "wb");
    assert(rec->spike_ids_file != NULL);
    
    rec->spike_times = calloc(size, sizeof(double));
    rec->spike_ids = calloc(size, sizeof(uint32_t));

    sprintf(output_path, "%s/%08dpos.bin", p->simulation.output_path, p->simulation.output_file_code);
    rec->position_file = fopen(output_path, "wb");
    assert(rec->position_file != NULL);
    
    if (p->simulation.stop_record_time != 0) {    
        rec->record_flag = true;

	sprintf(output_path, "%s/%08dvm.bin", p->simulation.output_path, p->simulation.output_file_code);
	rec->vm_file = fopen(output_path, "wb");
	assert(rec->vm_file != NULL);
	
	sprintf(output_path, "%s/%08dge.bin", p->simulation.output_path, p->simulation.output_file_code);
	rec->ge_file = fopen(output_path, "wb");
	assert(rec->ge_file != NULL);
	
	sprintf(output_path, "%s/%08dgi.bin", p->simulation.output_path, p->simulation.output_file_code);
	rec->gi_file = fopen(output_path, "wb");
	assert(rec->gi_file != NULL);

	rec->downsample_factor = p->simulation.record_downsample_factor;
	rec->start_bin = (int)floor(p->simulation.start_record_time/p->simulation.dt);
	rec->stop_bin = (int)floor(p->simulation.stop_record_time/p->simulation.dt);

	size_t rec_time_bins = (rec->stop_bin - rec->start_bin);
	rec->bin_size = ceil(rec_time_bins/p->simulation.record_downsample_factor);
	size = p->network.N * rec->bin_size;

	rec->vm = calloc(size, sizeof(float));
	rec->ge = calloc(size, sizeof(float));
        rec->gi = calloc(size, sizeof(float));
    }

    return NETSIM_NOERROR;
}

int initialize_partial_recorders(struct recorder *rec, struct parameters *p)
{
    /* set recording flag */
    rec->record_flag = false;
    rec->spike_count = 0;
    rec->buffer_pos = 0;
    
    rec->partial_time = p->simulation.partial_recording_time;
    rec->partial_spike_count = 0;
    rec->partial_file_count = 0;

    size_t size = p->network.N * 100 * rec->partial_time;
     
    rec->spike_times = calloc(size, sizeof(double));
    rec->spike_ids = calloc(size, sizeof(uint32_t));

    if (p->simulation.partial_recording_time != 0) {    
	rec->record_flag = true;

	rec->downsample_factor = p->simulation.record_downsample_factor;
	rec->start_bin = 0;
	rec->stop_bin = (int)floor(p->simulation.T/p->simulation.dt);
	
	size_t rec_time_bins = (rec->stop_bin - rec->start_bin)/rec->downsample_factor;
		
	rec->partial_bins = rec->partial_time/p->simulation.dt;
	rec->bin_size = ceil(rec->partial_bins/rec->downsample_factor);
	
	rec->partial_file_max = ceil(rec_time_bins/(double)rec->bin_size);
	rec->partial_bin_indices = calloc(rec->partial_file_max, sizeof(int));
	for(int i=0; i < rec->partial_file_max; i++) {
	    const int end = rec->start_bin + ((i+1) * rec->partial_bins);
	    const int last_bin = rec->stop_bin - rec->downsample_factor;
	    rec->partial_bin_indices[i] = (end > last_bin) ? last_bin : end;
	}
	
	size = p->network.N * rec->bin_size;
	
	rec->vm = calloc(size, sizeof(float));
	rec->ge = calloc(size, sizeof(float));
        rec->gi = calloc(size, sizeof(float));
    }

     return NETSIM_NOERROR;
}

int destory_recorders(struct recorder *rec, struct parameters *p)
{
#if !defined(LONG_SIMULATION)
    fclose(rec->spike_times_file);
    fclose(rec->spike_ids_file);
#endif

    /* free recorder*/
    free(rec->spike_times);
    free(rec->spike_ids);
    free(rec->position);

    if (p->simulation.stop_record_time != 0) {
#if !defined(LONG_SIMULATION)
	fclose(rec->vm_file);
	fclose(rec->ge_file);
	fclose(rec->gi_file);
#else
	free(rec->partial_bin_indices);
#endif
	free(rec->vm);
	free(rec->gi);
	free(rec->ge);
    }

    return NETSIM_NOERROR;
}

int record_spike(struct recorder *rec, const int N, const double T, const double t, const uint32_t ii)
{
#if !defined(LONG_SIMULATION)
    /* ensure spike count doesn't lead to overflow */
    assert(rec->spike_count < (N * T * 100));
	
    rec->spike_times[rec->spike_count] = t;
    rec->spike_ids[rec->spike_count] = ii;
    rec->spike_count++;
#else
    /* ensure spike count doesn't lead to overflow */
    assert(rec->partial_spike_count < (N * rec->partial_time * 100));
	
    rec->spike_times[rec->partial_spike_count] = t;
    rec->spike_ids[rec->partial_spike_count] = ii;
    rec->partial_spike_count++; 
    rec->spike_count++;  
#endif
    
    return NETSIM_NOERROR;
}

int record_state_variables(const struct simulation_parameters *sim,
			   struct recorder *rec, struct neuron *network,
			   const unsigned int N, const unsigned int tt)
{    
    if(rec->record_flag &&
       (tt % rec->downsample_factor) == 0 &&
       tt >= rec->start_bin &&
       tt < rec->stop_bin) {

	write_partial_variables(sim, rec, tt, N);	
	
	for(int ii = 0; ii < N; ii++) {
	    const int idx = (rec->buffer_pos * N) + ii;
	    rec->vm[idx] = network->v[ii];
	    rec->ge[idx] = network->ge[ii];
	    rec->gi[idx] = network->gi[ii];
	}

       	rec->buffer_pos++;
    }
    
    /* return */
    return NETSIM_NOERROR;
}

int record_static_variables(struct recorder *rec, struct neuron *network, const unsigned int N)
{
    rec->position = calloc(N * 3, sizeof(double));

    for(int ii = 0; ii < N; ii++) {
        rec->position[(ii * 3) + 0] = network->x_position[ii];
        rec->position[(ii * 3) + 1] = network->y_position[ii];
        rec->position[(ii * 3) + 2] = network->z_position[ii];
    }

    return NETSIM_NOERROR;
}

int write_partial_variables(const struct simulation_parameters *sim,
			    struct recorder *rec,
			    const unsigned int tt, const unsigned int N)
{
    
#if !defined(LONG_SIMULATION)
    return NETSIM_NOERROR;
#endif
    if(tt == rec->partial_bin_indices[rec->partial_file_count]) {
	char output_path[256];
    
	/* write out spikes and reset counter */
	sprintf(output_path, "%s/%08dspk_times_%04d.bin", sim->output_path,
		sim->output_file_code, rec->partial_file_count);
	rec->spike_times_file = fopen(output_path, "wb");
	assert(rec->spike_times_file != NULL);

	sprintf(output_path, "%s/%08dspk_ids_%04d.bin", sim->output_path,
		sim->output_file_code, rec->partial_file_count);
	rec->spike_ids_file = fopen(output_path, "wb");
	assert(rec->spike_ids_file != NULL);
	    
	fwrite(rec->spike_times, sizeof(double), rec->partial_spike_count, rec->spike_times_file);
	fwrite(rec->spike_ids, sizeof(uint32_t), rec->partial_spike_count, rec->spike_ids_file);

	fclose(rec->spike_times_file);
	fclose(rec->spike_ids_file);
	
	/* reset position in spike buffer*/
	rec->partial_spike_count = 0;

	if (rec->record_flag) {
	 /* write out state variables */
	 sprintf(output_path, "%s/%08dvm_%04d.bin", sim->output_path,
		 sim->output_file_code, rec->partial_file_count);
	 rec->vm_file = fopen(output_path, "wb");
	 assert(rec->vm_file != NULL);
	
	 sprintf(output_path, "%s/%08dge_%04d.bin", sim->output_path,
		 sim->output_file_code, rec->partial_file_count);
	 rec->ge_file = fopen(output_path, "wb");
	 assert(rec->ge_file != NULL);
	
	 sprintf(output_path, "%s/%08dgi_%04d.bin", sim->output_path,
		 sim->output_file_code, rec->partial_file_count);
	 rec->gi_file = fopen(output_path, "wb");
	 assert(rec->gi_file != NULL);

	 size_t size = N * rec->partial_bins/rec->downsample_factor;
	 fwrite(rec->vm, sizeof(float), size, rec->vm_file);
	 fwrite(rec->ge, sizeof(float), size, rec->ge_file);
	 fwrite(rec->gi, sizeof(float), size, rec->gi_file);

	 fclose(rec->vm_file);
	 fclose(rec->ge_file);
	 fclose(rec->gi_file);
	}

	/*reset buffer position*/
	rec->buffer_pos = 0;
	/* increment for next file chunk*/
	rec->partial_file_count++;
    
    }
     
    return NETSIM_NOERROR;
}

int write_output_files(const struct simulation_parameters *sim,
		       struct recorder *rec, const unsigned int N)
{
#if defined(LONG_SIMULATION)
    return NETSIM_NOERROR;
#endif
    
    /* write out static variables*/
    fwrite(rec->position, sizeof(double), N * 3, rec->position_file);
    fclose(rec->position_file);
    
    /* write spike file output */
    fwrite(rec->spike_times, sizeof(double), rec->spike_count, rec->spike_times_file);  
    fwrite(rec->spike_ids, sizeof(uint32_t), rec->spike_count, rec->spike_ids_file);
    
    /* write state variables */
    if (rec->record_flag) {
	fwrite(rec->vm, sizeof(float), N * rec->bin_size, rec->vm_file);
	fwrite(rec->ge, sizeof(float), N * rec->bin_size, rec->ge_file);
	fwrite(rec->gi, sizeof(float), N * rec->bin_size, rec->gi_file);
    }
    
    /* return */
    return NETSIM_NOERROR;
}
