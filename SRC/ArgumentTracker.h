//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#include <set>
template <class Enum>
class ArgumentTracker {
public:
  Enum current() const {
    int i = 0;
    while (static_cast<Enum>(i) != Enum::End
        && consumed.find(static_cast<Enum>(i)) != consumed.end()) {
      i++;
    }
    return static_cast<Enum>(i);
  }

  void consume(Enum argument) {
    consumed.insert(argument);
  }

  bool contains(Enum arg) {
    return consumed.find(arg) != consumed.end();
  }

  void increment() {
    consumed.insert(current());
  }

private:
  std::set<Enum>    consumed;
};

