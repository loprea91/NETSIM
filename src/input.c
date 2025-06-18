#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "input.h"

int parse_input_file(struct parameters *p) {
    /* init */
    FILE *file = NULL;
    char *line = NULL;
    char param[2048] = {0};
    char str[2048] = {0};
    double value = 0;
    int i = 0;
    int j = 0;
    size_t len = 0;

    /* open file */
    file = fopen(p->simulation.parameter_file_path, "r");
    if (file == NULL) {
        printf("Problem loading file\n\n");
        return NETSIM_ERROR;
    }

    /* read file line by line */
    while (getline(&line, &len, file) != EOF) {

        if ((*line != '#') && (*line != '$') && (*line != '%') && (*line != '\n') && (*line != ' '))
        /* line comment character: #, string char: $, linspace char: @ */
        {
            sscanf(line, "%s = %lf", param, &value);

            /* translation layer */
            if (strcmp(param, "N") == 0) {
                p->network.N = (int)round(value);
            } else if (strcmp(param, "K") == 0) {
                p->network.K = (int)round(value);
	    } else if (strcmp(param, "Nepop") == 0) {
                p->network.Nepop = value;
            } else if (strcmp(param, "Kin") == 0) {
                p->network.Kin = value;
            } else if (strcmp(param, "Kout") == 0) {
                p->network.Kout = value;
            } else if (strcmp(param, "L") == 0) {
                p->network.L = value;
            } else if (strcmp(param, "sigma_space") == 0) {
                p->network.sigma_space = value;
            } else if (strcmp(param, "poisson_rate") == 0) {
                p->network.poisson_rate = value;
            } else if (strcmp(param, "poisson_start_time") == 0) {
                p->network.poisson_start_time = value;
            } else if (strcmp(param, "poisson_stop_time") == 0) {
                p->network.poisson_stop_time = value;
            } else if (strcmp(param, "poisson_K_E") == 0) {
                p->network.poisson_K_E = value;
            } else if (strcmp(param, "poisson_K_I") == 0) {
                p->network.poisson_K_I = value;
            } else if (strcmp(param, "rewiring_probability") == 0) {
                p->network.rewiring_probability = value;
            }

            else if (strcmp(param, "T") == 0) {
                p->simulation.T = value;
            } else if (strcmp(param, "dt") == 0) {
                p->simulation.dt = value;
            } else if (strcmp(param, "vm_sigma") == 0) {
                p->simulation.vm_sigma = value;
            } else if (strcmp(param, "vm_mean") == 0) {
                p->simulation.vm_mean = value;
            } else if (strcmp(param, "ge_sigma") == 0) {
                p->simulation.ge_sigma = value;
            } else if (strcmp(param, "ge_mean") == 0) {
                p->simulation.ge_mean = value;
            } else if (strcmp(param, "gi_sigma") == 0) {
                p->simulation.gi_sigma = value;
            } else if (strcmp(param, "gi_mean") == 0) {
                p->simulation.gi_mean = value;
            } else if (strcmp(param, "start_record_time") == 0) {
                p->simulation.start_record_time = value;
            } else if (strcmp(param, "stop_record_time") == 0) {
                p->simulation.stop_record_time = value;
            } else if (strcmp(param, "record_downsample_factor") == 0) {
                p->simulation.record_downsample_factor = (int)round(value);
            } else if (strcmp(param, "connector") == 0) {
                p->simulation.connector = (int)round(value);
            } else if (strcmp(param, "initiator") == 0) {
                p->simulation.initiator = (int)round(value);
            } else if (strcmp(param, "output_file_code") == 0) {
                p->simulation.output_file_code = (int)round(value);
            } else if (strcmp(param, "connection_seed") == 0) {
                p->simulation.connection_seed = (int)round(value);
            } else if (strcmp(param, "conduction_speed_code") == 0) {
                p->simulation.conduction_speed_code = value;
            } else if (strcmp(param, "gamma_parameter_shape") == 0) {
                p->simulation.gamma_parameter_shape = value;
            } else if (strcmp(param, "gamma_parameter_scale") == 0) {
                p->simulation.gamma_parameter_scale = value;
            } else if (strcmp(param, "save_connectivity") == 0) {
                p->simulation.save_connectivity = value;
            } else if (strcmp(param, "save_positions") == 0) {
                p->simulation.save_positions = value;
            } else if (strcmp(param, "partial_recording_time") == 0) {
		p->simulation.partial_recording_time = value;
	    }

            else if (strcmp( param, "stim" ) == 0 ) {
                p->simulation.stim = value;  
            } else if (strcmp( param, "stim_id" ) == 0 ) {
                p->simulation.stim_id = value;  
            } else if (strcmp( param, "stim_time" ) == 0 ) {
                p->simulation.stim_time = value;  
            } else if (strcmp( param, "perturbation" ) == 0 ) {
                p->simulation.perturbation = value;  
            }

            else if (strcmp(param, "ge") == 0) {
                p->synapse.ge = value;
            } else if (strcmp(param, "gi") == 0) {
                p->synapse.gi = value;
            } else if (strcmp(param, "Ee") == 0) {
                p->synapse.Ee = value;
            } else if (strcmp(param, "Ei") == 0) {
                p->synapse.Ei = value;
            } else if (strcmp(param, "synapse_delay") == 0) {
                p->synapse.synapse_delay = value;
            } else if (strcmp(param, "min_conduction_speed") == 0) {
                p->synapse.min_conduction_speed = value;
            } else if (strcmp(param, "max_conduction_speed") == 0) {
                p->synapse.max_conduction_speed = value;
            } else if (strcmp(param, "taue") == 0) {
                p->synapse.taue = value;
            } else if (strcmp(param, "taui") == 0) {
                p->synapse.taui = value;
            } else if (strcmp(param, "min_uniform_delay") == 0) {
                p->synapse.min_uniform_delay = value;
            } else if (strcmp(param, "max_uniform_delay") == 0) {
                p->synapse.max_uniform_delay = value;
            } else if (strcmp(param, "p_release") == 0) {
                p->synapse.p_release = value;
            }

            else if (strcmp(param, "taum") == 0) {
                p->neuron.taum = value;
            } else if (strcmp(param, "vr") == 0) {
                p->neuron.vr = value;
            } else if (strcmp(param, "vreset") == 0) {
                p->neuron.vreset = value;
            } else if (strcmp(param, "vth") == 0) {
                p->neuron.vth = value;
            } else if (strcmp(param, "taur") == 0) {
                p->neuron.taur = value;
            } else if (strcmp(param, "El") == 0) {
                p->neuron.El = value;
            } else if (strcmp(param, "Ie") == 0) {
                p->neuron.Ie = value;
            } else if (strcmp(param, "Cm") == 0) {
                p->neuron.Cm = value;
            }

            else if (strcmp(param, "current_injection_count") == 0) {
                p->simulation.current_injection_count = value;
            }

            else {
                printf("Improper parameter read for %s.\n", param);
                return NETSIM_BAD_PARAMETER;
            }

        } else if (*line == '$') {

            /* sscanf */
            sscanf(line + 1, "%s = %s", param, str);
            printf("Reading %s: %s \n", param, str);

            /* translation layer */
            if (strcmp(param, "optional_connection_path") == 0) {
                strcpy(p->simulation.optional_connection_path, str);
            } else if (strcmp(param, "output_path") == 0) {
                if (p->simulation.output_path_override == 0) {
                    strcpy(p->simulation.output_path, str);
                }
            } else if (strcmp(param, "current_injection_path") == 0) {
                strcpy(p->simulation.current_injection_path, str);
            } else if (strcmp(param, "sparse_current_injection_path") == 0) {
                strcpy(p->simulation.sparse_current_injection_path, str);
            } else if (strcmp(param, "poisson_rate_path") == 0) {
                strcpy(p->simulation.poisson_rate_path, str);
            } else {
                printf("Improper parameter read (string).\n");
                return NETSIM_BAD_PARAMETER;
            }

        } else if (*line == '%') {
            /* retrieve param indices, and value*/
            sscanf(line + 1, "%2047[a-zA-Z_][%i, %i] = %lf", param, &i, &j, &value);

            /* assign values to matrices */
	}
    }

    /* cleanup */
    fclose(file);
    free(line);

    /* return */
    return NETSIM_NOERROR;
}

int load_array_from_file(double *array, int length, char *file_path) {
    FILE *file = fopen(file_path, "rb");
    if (file == NULL) {
        printf("Problem loading file\n\n");
        return NETSIM_BAD_PARAMETER;
    }

    /* read file */
    const int nbytes = fread(array, sizeof(double), length, file);
    if (nbytes != length) {
        printf("Problem reading file\n\n");
        fclose(file);
        return NETSIM_BAD_PARAMETER;
    }

    /* clean up */
    fclose(file);
    return NETSIM_NOERROR;
}

/* function: load file connections (binary) */
int file_connect(struct parameters *p, struct neuron *network) {
    /* init */
    const int K = p->network.K;
    FILE *file = fopen(p->simulation.optional_connection_path, "rb");

    /* loop over network */
    for (int ii = 0; ii < p->network.N; ii++) {
        if (fread((*(network + ii)).K, sizeof(uint32_t), K, file) != 0) {
            printf("Error in reading synaptic connections (ii=%d)\n", ii);
            return NETSIM_ERROR;
        }
    }

    /* clean up */
    fclose(file);
    return NETSIM_NOERROR;
}

int init_simulated_inputs(struct parameters *p) {
    if (strcmp(p->simulation.poisson_rate_path, "") != 0) {
        external_input_init(p);
    }
    if (strcmp(p->simulation.current_injection_path, "") != 0) {
        current_injection_init(p);
    }
    if (strcmp(p->simulation.sparse_current_injection_path, "") != 0) {
        sparse_current_injection_init(p);
    }

    return NETSIM_NOERROR;
}

int external_input_init(struct parameters *p) {
    int buffer_length = p->network.N * round(p->simulation.T / p->simulation.dt);
    p->simulation.poisson_variable_rate = calloc(buffer_length, sizeof(double));

    FILE *f = fopen(p->simulation.poisson_rate_path, "rb");
    assert(f != NULL);

    int r = fread(p->simulation.poisson_variable_rate, sizeof(double), buffer_length, f);
    assert(r == buffer_length);
    fclose(f);

    return NETSIM_NOERROR;
}

int current_injection_init(struct parameters *p) {
    int buffer_length =
        p->simulation.current_injection_count * round(p->simulation.T / p->simulation.dt);
    p->simulation.current_injection = calloc(buffer_length, sizeof(double));

    FILE *f = fopen(p->simulation.current_injection_path, "rb");
    assert(f != NULL);

    int r = fread(p->simulation.current_injection, sizeof(double), buffer_length, f);
    assert(r == buffer_length);

    fclose(f);

    return NETSIM_NOERROR;
}

int sparse_current_injection_init(struct parameters *p) {
    FILE *f = fopen(p->simulation.sparse_current_injection_path, "rb");
    assert(f != NULL);

    int n = 0;
    int r = fread(&n, sizeof(uint32_t), 1, f);
    assert(r == 1);
    
    p->simulation.sparse_current_entry_count = n;
    p->simulation.sparse_injection = calloc(n, sizeof(struct sparse_current_injection_entry));

    for (uint32_t i = 0; i < n; i++) {
        int r = fread(&p->simulation.sparse_injection[i].start_step, sizeof(uint32_t), 1, f);
        assert(r == 1);
        r = fread(&p->simulation.sparse_injection[i].stop_step, sizeof(uint32_t), 1, f);
        assert(r == 1);

        int count = 0;
        r = fread(&count, sizeof(uint32_t), 1, f);
        assert(r == 1);
        p->simulation.sparse_injection[i].index_count = count;
        p->simulation.sparse_injection[i].indices = calloc(count, sizeof(uint32_t));
        r = fread(p->simulation.sparse_injection[i].indices, sizeof(uint32_t), count, f);
        assert(r == count);
        r = fread(&p->simulation.sparse_injection[i].intensity, sizeof(double), 1, f);
        assert(r == 1);
    }
    fclose(f);

    return 0;
}
