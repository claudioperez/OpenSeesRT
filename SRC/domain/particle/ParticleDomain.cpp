#include <stdio.h>
#include <ParticleDomain.h>
  

ParticleDomain::ParticleDomain(int totnode, int maxfam)
  // Call the constructors for our member data
  : particles(totnode) // initialize the container of particles with size `totnode`
{

}


int ParticleDomain::hello(double x)
{
  // print the given double
  printf("hello, %lf\n", x);

  // return an integer
  return 32;
}

void ParticleDomain::print(int flag)
{
  printf("\"Particles\": [\n\t");
  for (Particle<3>& particle : particles) {
    particle.print(flag);
  }
  printf("]\n");
}
