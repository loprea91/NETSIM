#ifndef _INPUT_H_
#define _INPUT_H_

#include "simulation.h"

int parse_input_file(struct parameters *p);
int load_array_from_file(double *array, int length, char *file_path);

int init_simulated_inputs(struct parameters *p);
int external_input_init(struct parameters *p);
int current_injection_init(struct parameters *p);
int sparse_current_injection_init(struct parameters *p);

#endif
