#include <stdio.h>
#include <PeriDomain.h>
  

PeriDomain::PeriDomain(int totnode, int maxfam)
  // Call the constructors for our member data
  : particles(totnode) // initialize the container of particles with size `totnode`
{

}


int PeriDomain::hello(double x)
{
  // print the given double
  printf("hello, %lf\n", x);

  // return an integer
  return 32;
}

void PeriDomain::print(int flag)
{
  printf("  \"Particles\": [\n");
  for (PeriParticle<3>& particle : particles) {
    printf("    ");
    particle.print(flag);
  }
  printf("  ]\n");
}
