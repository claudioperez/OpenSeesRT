import re
import sys
import inspect
import textwrap
import importlib


from opensees import ast, obj


write = lambda *args, **kwds: print(*args, **kwds)


def get_type(a):
    pass

def write_grp(a):
    args = (write_grp(i) if isinstance(i, ast.Grp) else (i.name if i.name else "") for i in a.args)
    return "["+",".join(args)+"]"
    #write(f"<td>`{name}`</td><td>`typ`</td><td>{about}", end="")

def write_arg(a):
    write("<tr>",end="")
    typ = f"<code>{a.__class__.__name__}</code>"
    # about = a.about.replace('\n',' ').strip()
    name = a.name or ""
    default = " = "+str(a.default) if a.default else ""
    about = re.sub('[\s+]', ' ', a.about.replace('\n',' '))
    # write(f">| {a.name} | `{typ}` |  {about} |")
    # print(f">>>> {name} ")

    if isinstance(a, ast.Grp):
        if name == "":
            name = f"<code>{write_grp(a)}</code>"
        #elif a.type:
        #    typ = f"<code>{a.type.__name__}</code>"
        else:
            typ = f"<code>{write_grp(a)}</code>"
        write(f"<td>{name+default}</td><td>{typ}</td><td>{about}", end="")
        write("<table>")
        [write_arg(i) for i in a.args]
        write("</table>")

    elif isinstance(a, ast.Map):
        # typ = f"<code>{a.__class__.__name__}({a.key_type.__class__.__name__}: {a.val_type.__name__})</code>"
        typ = f"<code>{a.__name__}</code>"
        write(f"<td>{name+default}</td><td>{typ}</td><td>{about}", end="")
        write("<table>")
        write_arg(a.key_type)
        write_arg(a.val_type)
        write("</table>")

    else:
        if isinstance(a, ast.Ref):
            typ = f"<code>{a.__class__.__name__}({a.type.__name__})</code>"
        write(f"<td>{name+default}</td><td>{typ}</td><td>{about}", end="")

    write("</tr>")

def write_obj(v):
    s = str(inspect.signature(v)).replace('=None','')
    if len(s) > 65:
        #<span style="padding-left: 10%">
        s = s.replace(", ", ",<br>&emsp;&emsp;&emsp;")
    write(textwrap.dedent(f"""
    <!-- <blockquote> -->
    <span style="font-feature-settings: kern; color: var(--md-code-fg-color) !important; font-family: var(--md-code-font-family);">
        <span style="color:#900">{k}</span>{s}
    </span>
    """))
    try:
        write(v.about)
    except:
        write(v.__doc__ or "")
    if hasattr(v, "_img"):
        write(f"![](/figures/{v._img})")

    # write(">| | | |\n|--|--|--|")
    # write("<table>")
    write(textwrap.dedent("""
    <table>
    <colgroup>
      <col style="width: 10%" ><col style="width: 30%" ><col style="width: 60%" >
    </colgroup>
    <tbody>
    """))
    # write("<tbody>")
    try:
        pass
    except Exception as e:
        print(e, file=sys.stderr)
    if True:
        for a in v._args:
            # print(a, file=sys.stderr)
            write_arg(a)

    write(textwrap.dedent("""
    </tbody>
    </table>
    <!-- </blockquote> -->
    """))
    #write(":::\n")



def write_single(obj_name):
    lib = vars(getattr(lib, obj_name.split(".")[0]))
    write(textwrap.dedent(getattr(lib, attrs[0]).__doc__ or ""))


if __name__=="__main__":
    import opensees.lib


    libs = []
    attrs = []
    solo = False
    args = iter(sys.argv[1:])
    for arg in args:
        if arg[0] != "-":
            libs.append(arg)
        elif arg == "--attr":
            attrs.append(next(args))
        elif arg =="--single":
            solo = True

    if solo:
        lib = importlib.import_module(libs[0])
        write_single(lib, attrs[0])
        sys.exit()

    for lib_name in libs:
        lib = importlib.import_module(lib_name)
        head = (attrs or lib_name.split("."))[-1]
        write(textwrap.dedent(f"""
        ---
        title: {head}
        ...

        <style>
        h1 {{
            font-family: var(--md-code-font-family);
            color: var(--md-code-fg-color) !important;
            font-feature-settings: "kern";
        }}
        </style>

        # {head}

        """))
        if not attrs:
            objs = vars(lib)
            write(textwrap.dedent(lib.__doc__ or ""))

        else:
            objs = vars(getattr(lib, attrs[0]))
            write(textwrap.dedent(getattr(lib, attrs[0]).__doc__ or ""))

        write('<div style="width: 95%; padding-left: 5%">')

        for k,v in objs.items():
            if isinstance(v, type):
                if issubclass(v, obj.Component) and hasattr(v, "_args"):
                    write_obj(v)
            else:
                pass
        write("</div>")






