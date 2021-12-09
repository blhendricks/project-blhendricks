module tst
using BenchmarkTools

include("Event_Generator_symbols.jl")
import .Event_Gen_sym

#this describes the volume of the simulation
#vol = Dict("fiducial_rmin" => 0.0, "fiducial_rmax" => 5000.0, "fiducial_zmin" => -2700, "fiducial_zmax" => 0.0)

#num_list = [1, 10, 100, 1000]
#benchmarks_0 = benchmark_my_funciton(Event_gen.generate_eventlist_cylinder, num_list)

#plt0 = scatter(benchmarks_0.num_list, log10.(benchmarks_0.num_list./benchmarks_0.times_list), xscale=:log10, label=:none)
#		xlabel!(plt0, "Size of Array")
#		ylabel!(plt0, "log_10 (Evals/s)")
#		title!(plt0,"Benchmarks for my_function_0_args")
#		plt0

# The paramets in order of the following input are n_events, Emin, Emax,
# volume (Dict, found above), and interaction_type ("cc", "nc", or "ccnc")

#Event_Gen.generate_eventlist_cylinder(5, 1e19, 1e19, vol, "ccnc")
end
