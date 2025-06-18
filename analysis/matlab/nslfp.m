function lz = nslfp(params,v,ge,gi)

% make sure only the excitatory neurons are used
v = v( 1 : params.Ne, : );
ge = ge( 1 : params.Ne , : );
gi = gi( 1 : params.Ne , : );

tau = 0.006 / (params.dt * params.record_downsample_factor);
exc = circshift( ge .* (params.Ee - v), tau, 2 );
inh = gi.*(params.Ei - v);
lz = zscore( abs( exc - 1.65*inh ), [], 2 );

end