import opensees.tcl

rt = opensees.tcl.TclRuntime()

with open("./tests/models/meloland.M1.tcl","r") as f:
    rt.eval(f.read())

rt.eval("print -JSON -file /dev/stdout")

addr = rt._interp.tk.interpaddr()

print(f"{addr} ({hex(addr)})")

import libOpenSeesRT
builder = libOpenSeesRT.get_builder(addr)

print(builder.getSection("1"))

