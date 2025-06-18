#ifdef SPARSE_CURRENT_INJECTION

Ie = 0; /* default to zero */

for(int kk=0; kk < p->simulation.sparse_current_entry_count; kk++) {
	if(tt >= p->simulation.sparse_injection[kk].start_step &&
	   tt <  p->simulation.sparse_injection[kk].stop_step) {
		if(ii == p->simulation.sparse_injection[kk].indices[p->simulation.sparse_injection[kk].pos]) {
			Ie = p->simulation.sparse_injection[kk].intensity;
			p->simulation.sparse_injection[kk].pos++;
		}
		if(p->simulation.sparse_injection[kk].pos >=
		   p->simulation.sparse_injection[kk].index_count) {
			p->simulation.sparse_injection[kk].pos = 0;
		}
	}
}

#endif
