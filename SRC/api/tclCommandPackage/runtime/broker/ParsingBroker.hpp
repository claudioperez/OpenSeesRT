#include <string>

class FEM_BrokerAllClasses {

};


template <typename T, std::string name>
class ParsingBroker {

  ParsingBroker() {}

  T* readNewObject(std::string class_name, int argc, G3_Char** argv);

  T* makeNewObject(size_t hash) {return new T();}

  private:
    std::unordered_map<std::string, T*> m_object_map;
    std::string m_name = name;
};


