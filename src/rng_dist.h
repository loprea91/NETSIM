/********************************************************************************
 *                                                                              *
 *  FILE:     rng_dist.h                                                        *
 *  VERSION:  0.1                                                               *
 *  PROGRAM:  NETSIM                                                            *
 *                                                                              *
 *  PURPOSE:  Sampling from a different distributions 			        *
 *									        *
 *  DETAILS:  Returns a deviate distributed as a gamma distribution of integer  *
 * 	      order k, i.e., a waiting time to the kth event in a Poisson 	*
 * 	      process of unit mean, using UNIFORM_RAND as the source of         *
 * 	      uniform deviates. Adapted from Knuth, The Art of Computer         *
 * 	      Programming, vol2, 1997.					        *
 *
 *                                                                              *
 *  Copyright (C) 2016-2018 Lyle Muller                                         *
 *  http://cnl.salk.edu/~lmuller                                                *
 *                                                                              *
 * ---------------------------------------------------------------------------- *
 *                                                                              *
 *  DEVELOPMENT: Lyle Muller, Charlee Fletterman, Theo Desbordes                *
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

#ifndef RNG_DIST_H
#define RNG_DIST_H

#include "rng.h"

double gamma_sampling( int k, double scale, struct rng_state * rng );
double gamma_sampling_ab( double a, double b, struct rng_state * rng );
double beta_sampling( double a, double b, struct rng_state * rng  );
double _binomial_sampling( int t, double p, struct rng_state * rng );
double binomial_sampling_inversion( int t, double p, struct rng_state * rng );
#ifdef FAST_BINOMIAL
#define binomial_sampling( t, p, rng ) ( ((t) < 25) && ((t)>1e4) ? _binomial_sampling( t, p, rng ) : binomial_sampling_inversion( t, p, rng ) )
#else
#define binomial_sampling( t, p, rng ) _binomial_sampling( t, p, rng )
#endif

#endif /* RNG_DIST_H */
