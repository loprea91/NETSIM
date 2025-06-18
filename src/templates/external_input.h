#ifdef EXTERNAL_INPUT
if (tt >= poisson_start_time_steps &&
    tt < poisson_stop_time_steps) {
	const double poisson_spikes_E = binomial_sampling(poisson_K_E, dt * poisson_rate, rng) * we;
	buff_input_E[buffer_idx + ii] += poisson_spikes_E;

	const double poisson_spikes_I = binomial_sampling(poisson_K_I, dt * poisson_rate, rng) * wi;
	buff_input_I[buffer_idx + ii] += poisson_spikes_I;
}
#endif
