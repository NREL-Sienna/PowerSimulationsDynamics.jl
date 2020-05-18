using LITS
using PowerSystems
using Test
using NLsolve
using DiffEqBase
using Sundials
using InfrastructureSystems
using LinearAlgebra
using Logging

const IS = InfrastructureSystems
const PSY = PowerSystems

tests = readdir(dirname(@__FILE__))
tests = filter(
    f -> startswith(f, "test_") && endswith(f, ".jl") && f != basename(@__FILE__),
    tests,
)

const LOG_FILE = "lits-test.log"

LOG_LEVELS = Dict(
    "Debug" => Logging.Debug,
    "Info" => Logging.Info,
    "Warn" => Logging.Warn,
    "Error" => Logging.Error,
)

function get_logging_level(env_name::String, default)
    level = get(ENV, env_name, default)
    log_level = get(LOG_LEVELS, level, nothing)
    if isnothing(log_level)
        error("Invalid log level $level: Supported levels: $(values(LOG_LEVELS))")
    end

    return log_level
end

function run_tests()
    console_level = get_logging_level("SYS_CONSOLE_LOG_LEVEL", "Error")
    console_logger = ConsoleLogger(stderr, console_level)
    file_level = get_logging_level("SYS_LOG_LEVEL", "Info")

    #include("./data_tests/network_test_data.jl")
    include("./data_tests/dynamic_test_data.jl")
    include("./results/results_initial_conditions.jl")

    IS.open_file_logger(LOG_FILE, file_level) do file_logger
        multi_logger = IS.MultiLogger(
            [console_logger, file_logger],
            IS.LogEventTracker((Logging.Info, Logging.Warn, Logging.Error)),
        )
        global_logger(multi_logger)

        @testset "BasicTests" begin
            for test in tests
                print(splitext(test)[1], ": ")
                include(test)
                println()
            end
        end
        @info IS.report_log_summary(multi_logger)
    end
end

logger = global_logger()

try
    run_tests()
finally
    # Guarantee that the global logger is reset.
    global_logger(logger)
    nothing
end
