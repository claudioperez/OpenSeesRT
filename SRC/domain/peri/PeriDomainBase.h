#pragma once

// This class abstracts away the aspects of PeriDomain
// which are unrelated to the spatial dimension ndim.

class PeriDomainBase {
    public:
      PeriDomainBase(int totnode, int maxfam);

      int totnode, maxfam;

};