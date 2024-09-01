import opensees.tcl

def steel_cmd(name):
    fy  = 150e3
    E   = 100e3
    b   = 0.04
    R0  = 4
    cR1 = 0.9240
    cR2 = 0.1500
    return f"""
    model basic 1 1
    uniaxialMaterial Steel02 {name} {fy} {E} {b} {R0} {cR1} {cR2}
    """

def test_retrieve():

    rt = opensees.tcl.TclRuntime()
    rt.eval(steel_cmd("1"))
    rt.eval("print -json")

    # libOpenSeesRT must be imported by Python
    # AFTER if has been loaded by Tcl (this was done
    # when a TclRuntime() is created) so that Tcl stubs
    # are initialized. Otherwise there will be a segfault
    # when a python c-binding attempts to call a Tcl
    # C function. Users should never import libOpenSeesRT
    # themselves
    from opensees import OpenSeesPyRT

    addr = rt._interp.tk.interpaddr()

    builder = OpenSeesPyRT.get_builder(addr)

    print(builder.getUniaxialMaterial("1"))

def test_analyze():
    from opensees import uniaxial
    with uniaxial.Steel02(1, 150e3, 100e3, 0.04, 4.) as steel:
        steel.setTrialStrain(0.01)
        print(steel.getStress())

if __name__=="__main__":
    test_retrieve()
    test_analyze()


