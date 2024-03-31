#ifndef FrameTransform_h
#define FrameTransform_h
#include <CrdTransf.h>

template <int ndm, int nqv=3, int nn=2>
class FrameTransform : public CrdTransf {
  public:
  FrameTransform(int tag, int classTag) : CrdTransf(tag, classTag) {};
};
#endif
