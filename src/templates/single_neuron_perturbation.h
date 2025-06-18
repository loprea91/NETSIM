#ifdef SINGLE_NEURON_PERTURBATION

/* artificially spike one neuron in the network */
if (sim.stim == 1 && tt == sim.stim_time && ii == sim.stim_id) {
    printf("Setting neuron %d to threshold at time %f\n s", sim.stim_id, t);
    v[sim.stim_id] = vth;
}

#endif // SINGLE_NEURON_SPIKE