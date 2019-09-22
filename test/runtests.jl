using LITS
using PowerSystems
using Test

@test 1 == 1


tests = readdir(dirname( @__FILE__))
tests = filter(f -> startswith(f, "test_") &&
                               endswith(f, ".jl") &&
                               f != basename( @__FILE__),
                               tests)

for test in tests
    print(splitext(test)[1], ": ")
    include(test)
    println()
end

true
