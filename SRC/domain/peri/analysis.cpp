
#include <threads/thread_pool.hpp>
#include <mutex>

#include <cmath>
#include <stdio.h>
#include <vector>
#include <string>
#include <cstring>
#include <Logging.h>
#include <PeriDomain.h>
#include <PeriElement.h>
#include <PeriDomainBase.h>

#include <ElasticIsotropic.h>
#include <PeriParticle.h>
#include <NosbProj.h>



template <int ndim, int maxfam>
static int
PeriFormThreads(PeriDomain<ndim> &domain, std::vector<NosbProj<ndim, maxfam>>& nosb, OpenSees::thread_pool & threads)
{
    // OpenSees::thread_pool threads{6};
    std::mutex resp_mutex;

    // Form deformation gradients for trial
    threads.submit_loop<unsigned int>(0, nosb.size(), [&](int i) { 
        nosb[i].form_trial(); 
    }).wait();

    // Form force
    threads.submit_loop<unsigned int>(0, nosb.size(), [&](int i) {
			MatrixND<ndim, ndim> Q = nosb[i].sum_PKinv();

			for (int j = 0; j < nosb[i].numfam; j++)
			{
				const VectorND<ndim> T_j = nosb[i].bond_force(j, Q);

				const std::lock_guard<std::mutex> lock(resp_mutex);
				nosb[i].center->pforce += T_j;
				nosb[i].neigh[j]->pforce -= T_j;
			} 
    }).wait();
    return TCL_OK;
}