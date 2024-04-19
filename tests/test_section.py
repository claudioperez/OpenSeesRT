import sys
from opensees import tcl, section, patch

def test_section():
    print(tcl.ScriptBuilder(section.FiberSection(1, fibers=[patch.circ(extRad=0.3)])).script)

if __name__=="__main__":
    test_section()

