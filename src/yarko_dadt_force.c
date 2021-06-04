/**
 * @file    radiation_forces.c
 * @brief   Add radiation forces
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Yarkovsky Forces$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 R. Melikyan
 * Implementation Paper    TBD
 * Based on                `Granvik et al. 2017 <https://www.aanda.org/articles/aa/abs/2017/02/aa29252-16/aa29252-16.html>`_.
 * C Example               TBD
 * Python Example          TBD
 * ======================= ===============================================
 * 
 * This applies a semi-major axis drift to particles in the simulation.
 * This is the simplist implimentation of a transverse diurnal Yarkovsky effect.
 * Only particles whose `obliquity` parameter is set will feel the drift.  
 * 
 * **Effect Parameters**
 * 
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * dadt_1km (double)            Yes         Reference dadt heliocentric drift for a 1km diameter object in AU/Myr
 * ============================ =========== ==================================================================
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * primary (int)                No          Flag central body. Assumes index 0 
 * obliquity (double)           Yes         Object Spin angle in radians. Yarkovsky transverse drift will not be computed for objects without this parameter.
 * ============================ =========== ==================================================================
 * 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "reboundx.h"

static void rebx_calculate_yarko_dadt_force(struct rebx_extras* const rebx, struct reb_simulation* const sim, const double dadt_1km, const double star_i, struct reb_particle* const particles, const int N){
    const struct reb_particle star = particles[star_i];
    const double mu = sim->G*star.m;
    const double au2m = 149597870700;
    double sec2year = 31536000.0;


    for (int i=0;i<N;i++){
        
        const double* obliquity = rebx_get_param(rebx, particles[i].ap, "obliquity");
        if(obliquity == NULL) continue; // only particles with obliquity set feel yarko drift

        const double D = 2 * particles[i].r * 1e-3 // Particle Diameter in km
        
        const struct reb_particle p = particles[i];
        const double dx = p.x - star.x; 
        const double dy = p.y - star.y;
        const double dz = p.z - star.z;
        const double dr = sqrt(dx*dx + dy*dy + dz*dz); // distance to star
        
        const double dvx = p.vx - star.vx;
        const double dvy = p.vy - star.vy;
        const double dvz = p.vz - star.vz;
        const double v2   = dvx*dvx + dvy*dvy + dvz*dvz;

        const double energy     = 0.5*v2 - mu/dr;
        const double a          = -0.5*mu/energy;
        const double auPmyr2mPs = au2m/sec2year/1e6; // conversion factor for au/Myr to m/s
        const double dadt       = dadt_1km/D*cos(obliquity); // dadt for object of diameter D and obliquity
        const double k          = 0.5*dadt*auPmyr2mPs*mu/(a*a); 

        particles[i].ax += k * dvx/v2;
        particles[i].ay += k * dvy/v2;
        particles[i].az += k * dvz/v2;
	}
}

void rebx_yarko_dadt_force(struct reb_simulation* const sim, struct rebx_force* const radiation_forces, struct reb_particle* const particles, const int N){
    struct rebx_extras* const rebx = sim->extras;
    double* dadt_1km = rebx_get_param(rebx, radiation_forces->ap, "dadt_1km");
    if (dadt_1km == NULL){
        reb_error(sim, "Need to set a reference dadt drift for a 1km diameter object.  See examples in documentation.\n");
    }
    
    int star_found=0;
    for (int i=0; i<N; i++){
        if (rebx_get_param(rebx, particles[i].ap, "primary") != NULL){
            star_found = 1;
            rebx_calculate_yarko_dadt_force(rebx, sim, *dadt_1km, i, particles, N);
        }
    }
    if (!star_found){
        rebx_calculate_yarko_dadt_force(rebx, sim, *dadt_1km, 0, particles, N);    // default source to index 0 if "primary" not found on any particle
    }
}
