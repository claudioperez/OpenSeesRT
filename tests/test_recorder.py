import opensees
import quakeio.elcentro

def test_node():
    k = 10.
    m = 1.0
    model = opensees.model(ndf = 1,
        nodes = {1: [0], 2: [0]},
        elems = [[
                opensees.element.ZeroLength(
                    materials = [opensees.uniaxial.Elastic(E=k)]
                ),
                [1, 2]
            ]
        ],
        zeros = [
            [[1], [1]]
        ]
    )

    model.mass(2, [m])

    patterns = {
        1: opensees.pattern.UniformExcitation(dof=1, accel=quakeio.elcentro.accel)
    }

    filename = "/dev/stdout"
# filename = "a.txt"

    recorders = [
        opensees.recorder.Node(filename, recorder="disp", format="txt", nodes=[2], dofs=range(1,2), time=True)
    ]

    print(opensees.tcl.dumps(recorders[0]))

    analysis = opensees.analysis.DirectIntegrationAnalysis(model, patterns, recorders=recorders)

    for i in range(10):
        analysis.analyze(1, time_step=0.02)


# def test_fiber_beam():
#     opensees.recorder.Element("out.txt", {"section": [2], "fiber": [3]})
#     opensees.recorder.Element("out.txt", "")
test_node()


