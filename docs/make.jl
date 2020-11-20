using Documenter, PowerSystems, DocStringExtensions, PowerSimulationsDynamics

makedocs(
    modules = [PowerSimulationsDynamics],
    format = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename = "PowerSimulationsDynamics.jl",
    pages = Any[ # Compat: `Any` for 0.4 compat
        "Welcome Page" => "index.md",
        # "User Guide" => "man/guide.md",
        "Tutorials" => "tutorials_page.md",
        "Models" => Any[
            "Network" => "Models/network.md",
            "Reference Frames" => "Models/srf.md",
            "Generator" => "Models/gens.md",
            "Inverter" => "Models/inverters.md",
            "Small Signal" => "Models/small.md",
        ],
        "Public API Reference" => "api/public.md",
        "Internal API Reference" => "api/internal.md"
    ],
)

deploydocs(
    repo = "github.com/NREL-SIIP/PowerSimulationsDynamics.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
    push_preview=true,
    versions = ["stable" => "v^", "v#.#"]
)
