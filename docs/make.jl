using Documenter, PowerSimulationsDynamics

makedocs(
    modules = [PowerSimulationsDynamics],
    format = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename = "PowerSimulationsDynamics.jl",
    pages = Any[ # Compat: `Any` for 0.4 compat
        "Home" => "index.md",
        # "User Guide" => "man/guide.md",
        "Tutorials" => Any[
            "Tutorial 0: Data Creation" => "Examples/example_data.md",
            "Tutorial 1: OMIB" => "Examples/example_OMIB.md",
            "Tutorial 2: Dynamic Lines" => "Examples/example_lines.md",
        ],
        "Models" => Any[
            "Network" => "Models/network.md",
            "Reference Frames" => "Models/srf.md",
            "Generator" => "Models/gens.md",
            "Inverter" => "Models/inverters.md",
            "Small Signal" => "Models/small.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/NREL-SIIP/PowerSimulations.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
    versions = ["stable" => "v^", "v#.#"],
)
