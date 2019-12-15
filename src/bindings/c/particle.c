
#include "particle.h"

void initialise_particles_c(int *fileid);
void controller_particlesolve_c(int *fileid);


void initialise_particles(int *fileid)
{
  initialise_particles_c(fileid);
}

void controller_particlesolve(int *fileid)
{
  controller_particlesolve_c(fileid);
}
