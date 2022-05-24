using Documenter, PowerSystems, DocStringExtensions, PowerSimulationsDynamics

makedocs(
    modules = [PowerSimulationsDynamics],
    format = Documenter.HTML(mathengine = Documenter.MathJax()),
    sitename = "PowerSimulationsDynamics.jl",
    pages = Any[
        "Welcome Page" => "index.md",
        "Quick Start Guide" => "quick_start_guide.md",
        "Simulation Execution" => "execute.md",
        "Tutorials" => "tutorials_page.md",
        "Models" => "models.md",
        "Initialization" => "initialization.md",
        "Small Signal" => "small.md",
        "Reference Frames" => "reference_frames.md",
        "Perturbations" => "perturbations.md",
        "Industrial Renewable Models" => "generic.md",
        "Generator Component Library" => Any[
            "Machine" => "component_models/machines.md",
            "Shaft" => "component_models/shafts.md",
            "AVR" => "component_models/avr.md",
            "PSS" => "component_models/pss.md",
            "Turbine and Governor" => "component_models/turbine_gov.md",
        ],
        "CIG Component Library" => Any[
            "Converter" => "component_models/converter.md",
            "DC Sources" => "component_models/dc_source.md",
            "Filter" => "component_models/filters.md",
            "Frequency Estimators" => "component_models/freq_esti.md",
            "Inner Control" => "component_models/inner_control.md",
            "Outer Control" => "component_models/outer_control.md",
        ],
        "Branch Models" => Any["Network" => "component_models/network.md",],
        "Load Models" => Any["Load Models" => "component_models/loads.md",],
        "Code Base Developer Guide" => Any["Developer Guide" => "code_base_developer_guide/developer.md",],
        "Public API Reference" => "api/public.md",
        "Internal API Reference" => "api/internal.md",
    ],
)

deploydocs(
    repo = "github.com/NREL-SIIP/PowerSimulationsDynamics.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"],
)
