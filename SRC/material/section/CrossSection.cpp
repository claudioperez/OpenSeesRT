
constexpr unsigned long log2(unsigned long n)
{
  return ( (n<2) ? 1 : 1+log2(n/2));
}


enum FrameType : unsigned long {
  N     = 1<< 0, // 0b00000001
  Vy    = 1<< 1, // 0b00000010
  Vz    = 1<< 2, // 0b00000100
  T     = 1<< 3, // 0b00000000
  My    = 1<< 4, // 0b00000000
  Mz    = 1<< 5, // 0b00000000
  R     = 1<< 6, // 0b00000000
  Q     = 1<< 7, // 0b00000000
  B     = 1<< 8, // 0b00000000
  W     = 1<< 9, // 0b00000000
  End   = 1<<10
};

class FrameSection {
};


#if 0
enum class ShellType {
  N, Vy, Vz,
  T, My, Mz
};

class TrussSection {
};

class ShellSection {
};

#endif
