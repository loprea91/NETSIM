/********************************************************************************
 *                                                                              *
 *  FILE:     rng_dist.c                                                        *
 *  VERSION:  0.1                                                               *
 *  PROGRAM:  NETSIM                                                            *
 *                                                                              *
 *  PURPOSE:  Sampling from a PRNG distributions    			        *
 *										*
 *  DETAILS:  Returns a deviate distributed as a gamma distribution of integer  *
 * 	      order k, i.e., a waiting time to the kth event in a Poisson 	*
 * 	      process of unit mean, using UNIFORM_RAND as the source of 	*
 * 	      uniform deviates. Adapted from Knuth, The Art of Computer         *
 * 	      Programming, Volume 2, 1997.				        *
 *            The a,b paramatarization from Marsaglia and Tsang -               *
 *            "A Simple Method for Generating Gamma Variables" (2000)           *
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

#include <math.h>
#include <stdio.h>

#include "rng_dist.h"

/* function: gamma_sampling */
double gamma_sampling(int k, double scale, struct rng_state *rng) {
    int ii;
    double x;

    if (k < 1) {
        // printf( "Error in routine gamma sampling: unsupported parameter k (shape of the gamma
        // distribution)\n" );
    }

    /* Use direct method, adding waiting times. */
    x = 1.0;
    for (ii = 1; ii <= k; ii++) {
        x *= rng_dbl64(rng);
    }
    x = -log(x);

    return x * scale;
}

/* function: gamma_sampling_ab */
double gamma_sampling_ab(double a, double b, struct rng_state *rng) {
    double d, c, x, v, U;

    d = a - (1 / 3);
    c = 1 / sqrt(9 * d);

    while (1) {
        do {
            x = rng_gauss(rng);
            v = pow(1 + c * x, 3);
        } while (v <= 0);
        U = rng_dbl64(rng);
        if (U < 1 - 0.0331 * pow(x, 4)) {
            break;
        }
        if (log(U) < 0.5 * x * x + d * (1 - v + log(v))) {
            break;
        }
    }
    return (d * v) / b;
}

/* function: beta_sampling */
double beta_sampling(double a, double b, struct rng_state *rng) {
    double X1 = gamma_sampling_ab(a, 1, rng);
    double X2 = gamma_sampling_ab(b, 1, rng);
    if (a < 1.0 && b < 1.0) {
        while (X1 + X2 <= 1) {
            X1 = exp(log(gamma_sampling_ab(a, 1, rng)) / a);
            X1 = exp(log(gamma_sampling_ab(b, 1, rng)) / b);
        }
    }

    return X1 / (X1 + X2);
}

/* function: binomial_sampling */
double _binomial_sampling(int t, double p, struct rng_state *rng) {
    int a = floor(t / 2) + 1;
    int b = t - a + 1;

    double X = beta_sampling(a, b, rng);
    if (X >= p) {
        t = a - 1;
        p = p / X;
    } else {
        t = a + b - 1;
        p = (p - X) / (1 - X);
    }

    int count = 0;
    for (int i = 0; i < t; i++) {
        if (rng_dbl64(rng) < p) {
            count++;
        }
    }
    return count;
}

/* function: binomial_sampling_inversion */
double binomial_sampling_inversion(int n, double p, struct rng_state *rng) {
    int k = 0;
    double q = 1.0 - p; // number of successes
    double s = p / q; // ratio of successes to failures
    double a = (n + 1) * s; // number of expected successes
    double r = exp(n * log(q)); // inverse CDF @ p
    double u = rng_open_dbl(rng) ; // uniform random number on [0, 1)
    while (u > r) { // while u > inverse CDF @ p
        u -= r; // u is now a uniform random number on [0, r]
        k++; // increment number of successes
        r *= (a / k) - s; // update inverse CDF @ p to account for new success
    } // end while

    return k;
}
