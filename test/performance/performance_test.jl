precompile = @timed using PowerSimulationsDynamics

using PowerSimulationsDynamics
const PSID = PowerSimulationsDynamics
using Sundials
using PowerSystems
const PSY = PowerSystems
using OrdinaryDiffEq
using Logging

open("precompile_time.txt", "a") do io
    write(io, "| $(ARGS[1]) | $(precompile.time) |\n")
end

raw_file_dir = "test/data_tests/240busWECC_2018_PSS32_fixed_shunts.raw"
dyr_file = "test/data_tests/240busWECC_2018_PSS.dyr"
sys = System(
    raw_file_dir,
    dyr_file;
    bus_name_formatter = x -> string(x["name"]) * "-" * string(x["index"]),
)

# First runs

for m in [ResidualModel, MassMatrixModel]
    Simulation(
        m,
        sys,
        pwd(),
        (0.0, 20.0), #time span
        BranchTrip(1.0, Line, "CORONADO    -1101-PALOVRDE    -1401-i_10");
        console_level = Logging.Error,
    )
end

try
    sim_ida, time_build_ida, _, _ = @timed Simulation(
        ResidualModel,
        sys,
        pwd(),
        (0.0, 20.0), #time span
        BranchTrip(1.0, Line, "CORONADO    -1101-PALOVRDE    -1401-i_10");
        console_level = Logging.Error,
    )
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build ResidualModel | $(time_build_ida) |\n")
    end
catch e
    @error exception = (e, catch_backtrace())
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build ResidualModel | FAILED TO BUILD |\n")
    end
end

try
    sim_rodas, time_build_rodas, _, _ = @timed Simulation(
        MassMatrixModel,
        sys, #system
        pwd(),
        (0.0, 20.0), #time span
        BranchTrip(1.0, Line, "CORONADO    -1101-PALOVRDE    -1401-i_10");
        console_level = Logging.Error,
    ) #Type of Fault
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build MassMatrixModel | $(time_build_rodas) |\n")
    end
catch e
    @error exception = (e, catch_backtrace())
    open("execute_time.txt", "a") do io
        write(io, "| $(ARGS[1])-Build MassMatrixModel | FAILED TO BUILD |\n")
    end
end

#=
try
    status = execute!(sim_rodas, Rodas4(), dtmax = 0.01)
    res_rodas = read_results(sim_rodas)
catch
    time = "EXECUTE FAILED"
finally
    open("execute_time.txt", "a") do io
    write(io, "| $(ARGS[1]) | 999 |\n")
    end
end

try
    status = execute!(sim_ida, IDA(), dtmax = 0.01)
    res_ida = read_results(sim_ida)
catch
    time = "EXECUTE FAILED"
finally
    open("execute_time.txt", "a") do io
    write(io, "| $(ARGS[1]) | 999 |\n")
    end
end
=#
