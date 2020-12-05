#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "V_header.c"
#include "mtwister.h"


// --- data structures ---

typedef struct {
    double sqr_coeff;
    double inter_coeff;
    double V_tilde_coeff;
} action_coefficients_t;

typedef struct {
    unsigned long seed;
    unsigned int N;
    unsigned int local_steps;
    unsigned int measuration_to_take;
    unsigned int measure_every;
    double step_size;
} markov_setup_t;

typedef double* state_t;
typedef unsigned int* geometry_t; 

// --- parameters interpretes ---

// invocation_string = "{simulator_file_name} {seed} {d} {lambda} {eta} {N} {local_steps} {measuration_to_take} {measure_every} {step_size}"

action_coefficients_t get_coefficients(char **argv){
    double d = strtod(argv[2], NULL);
    double lambda = strtod(argv[3], NULL);
    double eta = strtod(argv[4], NULL);
    
    action_coefficients_t coeffs = {
        sqrt(d) * (1. / eta + eta * lambda / 2.),
        - sqrt(d) / eta,
        sqrt(d) * eta
    };
    return coeffs;
}

markov_setup_t get_markov_setup(char **argv){
    markov_setup_t setup = {
        strtoul(argv[1], NULL, 10),
        (unsigned int)strtoul(argv[5], NULL, 10),
        (unsigned int)strtoul(argv[6], NULL, 10),
        (unsigned int)strtoul(argv[7], NULL, 10),
        (unsigned int)strtoul(argv[8], NULL, 10),
        strtod(argv[9], NULL)
    };
    return setup; 
}

// --- main loop ---

// generate start state, using cold-start
state_t gen_start_state(unsigned int N){
    state_t new_state = (state_t) malloc(N * sizeof(double));
    for(int n = 0; n < N; n++)
        new_state[n] = 0;
    return new_state;
}

// generate geometry array containing the loop around
geometry_t gen_geometry(unsigned int N){
    geometry_t geometry = (unsigned int*)malloc(2*N*sizeof(unsigned int));
    geometry[0] = N - 1;
    geometry[N] = 1;
    for(unsigned int n = 1; n < N - 1; n++){
        geometry[n] = n - 1;
        geometry[N+n] = n + 1;
    }
    geometry[N - 1] = N - 2;
    geometry[2*N - 1] = 0;
    return geometry;
}

void run_simulation(action_coefficients_t coeffs, markov_setup_t setup){
    state_t state = gen_start_state(setup.N);  // generate new state
    geometry_t geometry = gen_geometry(setup.N);  // generate geometry
    MTRand rand_state = seedRand(setup.seed);  // initialize random number generator
    for(int measure = 0; measure < setup.measuration_to_take; measure++){
        unsigned int accepted = 0;  // count accepted ratio
        for(int passes = 0; passes < setup.measure_every; passes++)
            for(int n = 0; n < setup.N; n++)
                for(int local_step = 0; local_step < setup.local_steps; local_step++){
                    // do a local step on y_n
                    double delta = setup.step_size * (2*genRand(&rand_state) - 1);  // choose a point in [-step_size, step_size)
                    double delta_S = 
                        coeffs.sqr_coeff * (2*state[n]*delta + delta*delta) + //quadratic term
                        coeffs.inter_coeff * (state[geometry[n]] + state[geometry[setup.N + n]]) * delta + // interaction term
                        coeffs.V_tilde_coeff * (V_tilde(state[n] + delta) - V_tilde(state[n]));  // potential change
                    if(
                        delta_S < 0 ||  // delta_S < 0 is going to a more probable state, always accepted 
                        genRand(&rand_state) < exp(-delta_S)  // accepted with probability exp(-delta_S)
                    ){
                        state[n] += delta;  // changing state
                        accepted ++;
                    }
                }
        // print accepted ratio
        printf("%f", accepted / ((float) setup.measure_every*setup.N*setup.local_steps));
        // print out the state
        for(int n = 0; n < setup.N; n++)
            printf(" %f", state[n]);
        printf("\n");
    }
}

// --- main ---    

int main(int argc, char **argv){
    if(argc != 10){
        fprintf(stderr, "Wrong argument number");
        return 1;
    } 
    action_coefficients_t coeffs = get_coefficients(argv);
    markov_setup_t setup = get_markov_setup(argv);
    run_simulation(coeffs, setup);
    return 0;
}