#pragma once

// This class abstracts away the aspects of PeriDomain
// which are unrelated to the spatial dimension ndim.

class PeriDomainBase {
    public:
      PeriDomainBase(int totnode, int maxfam);

      // virtual destructor
      // default: explicitly instructs the compiler to generate
      // the default implementation of the destructor
      virtual ~PeriDomainBase() = default;
      
      // pure virtual function
      // a function that has no implementation in the base class
      // and must be overridden in the derived class
      virtual int getNDM() const = 0; // function to get the number of dimensions


      int totnode, maxfam;
      char plane_type;
      double delta = 0.0, space = 0.0;

};