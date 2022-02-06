
#ifdef _PARALLEL_PROCESSING
int
partitionModel(int eleTag)
{
  if (OPS_PARTITIONED == true)
    return 0;

  int result = 0;

  if (OPS_theChannels != 0)
    delete[] OPS_theChannels;

  OPS_theChannels = new Channel *[OPS_NUM_SUBDOMAINS];

  // create some subdomains
  for (int i = 1; i <= OPS_NUM_SUBDOMAINS; i++) {
    if (i != OPS_MAIN_DOMAIN_PARTITION_ID) {
      ShadowSubdomain *theSubdomain =
          new ShadowSubdomain(i, *OPS_MACHINE, *OPS_OBJECT_BROKER);
      theDomain.addSubdomain(theSubdomain);
      OPS_theChannels[i - 1] = theSubdomain->getChannelPtr();
    }
  }

  // create a partitioner & partition the domain
  if (OPS_DOMAIN_PARTITIONER == 0) {
    //      OPS_BALANCER = new ShedHeaviest();
    OPS_GRAPH_PARTITIONER = new Metis;
    // OPS_DOMAIN_PARTITIONER = new DomainPartitioner(*OPS_GRAPH_PARTITIONER,
    // *OPS_BALANCER);
    OPS_DOMAIN_PARTITIONER = new DomainPartitioner(*OPS_GRAPH_PARTITIONER);
    theDomain.setPartitioner(OPS_DOMAIN_PARTITIONER);
  }

  // opserr << "commands.cpp - partition numPartitions: " << OPS_NUM_SUBDOMAINS
  // << endln;

  result = theDomain.partition(OPS_NUM_SUBDOMAINS, OPS_USING_MAIN_DOMAIN,
                               OPS_MAIN_DOMAIN_PARTITION_ID, eleTag);

  if (result < 0)
    return result;

  OPS_PARTITIONED = true;

  DomainDecompositionAnalysis *theSubAnalysis;
  SubdomainIter &theSubdomains = theDomain.getSubdomains();
  Subdomain *theSub = 0;

  // create the appropriate domain decomposition analysis
  while ((theSub = theSubdomains()) != 0) {
    if (the_static_analysis != 0) {
      theSubAnalysis = new StaticDomainDecompositionAnalysis(
          *theSub, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
          *theSOE, *theStaticIntegrator, theTest, false);

    } else {
      theSubAnalysis = new TransientDomainDecompositionAnalysis(
          *theSub, *theHandler, *theNumberer, *the_analysis_model, *theAlgorithm,
          *theSOE, *theTransientIntegrator, theTest, false);
    }
    theSub->setDomainDecompAnalysis(*theSubAnalysis);
  }

  return result;
}

#endif
