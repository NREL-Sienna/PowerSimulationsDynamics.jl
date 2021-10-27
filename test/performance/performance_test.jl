res = @timed using PowerSimulationsDynamics

open("precompile_time.txt", "a") do io
    write(io, "| $(ARGS[1]) | $(res.time) |\n")
end

open("execute_time.txt", "a") do io
    write(io, "| $(ARGS[1]) | $(res.time) |\n")
end
