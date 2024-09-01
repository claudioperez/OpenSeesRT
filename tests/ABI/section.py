import opensees.tcl

# opensees.loadtcl(filename)

rt = opensees.tcl.TclRuntime()

#with open(".cache/meloland.M1.tcl","r") as f:
with open("./docs/gallery/sections/example_9/Ex9b.build.WSection2D.tcl","r") as f:
    rt.eval(f.read())

rt.eval("print -JSON -file /dev/stdout")

addr = rt._interp.tk.interpaddr()

print(f"{addr} ({hex(addr)})")

import libOpenSeesRT
builder = libOpenSeesRT.get_builder(addr)

print(builder.getSection("1"))

