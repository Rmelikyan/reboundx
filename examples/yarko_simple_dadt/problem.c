/**
 * Radiation forces on circumplanetary dust
 * 
 * This example shows how to integrate circumplanetary
 * dust particles under the action of radiation forces using IAS15.
 * We use Saturn's Phoebe ring as an example, a distant ring of debris, 
 * The output is custom, outputting the semi-major axis of 
 * every dust particle relative to the planet. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void heartbeat(struct reb_simulation* r);

const double sec2year = 31536000.0;
const double au2m = 149597870700;
const double tmax = 1; // max time in years

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    // Setup constants 
    sim->integrator     = REB_INTEGRATOR_WHFAST;
    sim->G              = 6.674e-11;    // Use SI units
    sim->dt             = 1e6;          // Initial timestep in sec
    sim->N_active       = 1;            // Only the sun affects other particles gravitationally
    sim->heartbeat      = heartbeat;
    
    // sun
    struct reb_particle sun = {0};
    sun.m  = 1.99e30;                   // mass of Sun in kg
    reb_add(sim, sun);
    
    // Add REBOUNDx
    struct rebx_extras* rebx = rebx_attach(sim); 
    struct rebx_force* rad = rebx_load_force(rebx, "yarko_dadt_force");
    double dadt_1km = 2e-4;             // Reference drift rate in au/Myr for 1km diameter body
    rebx_set_param_double(rebx, &rad->ap, "dadt_1km", dadt_1km);
    
    // Will assume particles[0] is the primary body by default. You can also add a flag to a particle explicitly
    rebx_set_param_int(rebx, &sim->particles[0].ap, "primary", 1); 

    /* Test particles
     Here we imagine particles launched from Saturn's irregular Satellite Phoebe.
     Such grains will inherit the moon's orbital elements (e.g. Tamayo et al. 2011) 
        
     In order for a particle to feel the yarkovsky drift force, we have to set their obliquity parameter, 
     the angle of their spin axis with respect to their orbiting plane. */
        
    double a_test = 2.5*au2m;            // semimajor axis of satellite Phoebe, in m
    double e_test = 0.1;               // eccentricity of Phoebe
    double inc_test = 0.;   // inclination of Phoebe to Saturn's orbital plane
    double Omega_test = 0.;             // longitude of ascending node
    double omega_test = 0.;             // argument of pericenter
    double f_test = 0.;                 // true anomaly

    // We first set up the orbit and add the particles
    double m_test = 0.;                 // treat dust particles as massless
    struct reb_particle p = reb_tools_orbit_to_particle(sim->G, sim->particles[0], m_test, a_test, e_test, inc_test, Omega_test, omega_test, f_test);

    for (int i=0;i<5;i++){
        p.r = i; // size determines drift rate
        reb_add(sim, p); 
        // assign particles obliquity of 0 or 180 degrees to demonstrate drift directions
        rebx_set_param_double(rebx, &sim->particles[-1].ap, "obliquity", i*M_PI);
    }
    
    reb_move_to_com(sim);

    reb_integrate(sim, tmax*sec2year);
    rebx_free(rebx);                /* free memory allocated by REBOUNDx */
}

void heartbeat(struct reb_simulation* sim){
    if(reb_output_check(sim, 1e-1*sec2year)){
        const struct reb_particle* particles = sim->particles;
        const double t = sim->t;
        struct reb_orbit orbit = reb_tools_particle_to_orbit(sim->G, particles[1], particles[0]); /* calculate orbit of particles[2] around Saturn */
        double delta_a1 = orbit.a/au2m - 2.5;
        orbit = reb_tools_particle_to_orbit(sim->G, particles[2], particles[0]); 
        double delta_a2 = orbit.a/au2m - 2.5;
        orbit = reb_tools_particle_to_orbit(sim->G, particles[3], particles[0]); 
        double delta_a3 = orbit.a/au2m - 2.5;
        orbit = reb_tools_particle_to_orbit(sim->G, particles[4], particles[0]); 
        double delta_a4 = orbit.a/au2m - 2.5;
        orbit = reb_tools_particle_to_orbit(sim->G, particles[5], particles[0]); 
        double delta_a5 = orbit.a/au2m - 2.5;
        printf("%e\tda1: %f\tda2: %f\tda3: %f\tda4: %f\tda5: %f\n", t/sec2year, delta_a1, delta_a2, delta_a3, delta_a4, delta_a5);
    }
}
