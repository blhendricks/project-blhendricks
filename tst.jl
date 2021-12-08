include("Event_Generator.jl")
import .Event_Gen

#this describes the volume of the simulation
vol = Dict("fiducial_rmin" => 0, "fiducial_rmax" => 4, "fiducial_zmin" => -2.7, "fiducial_zmax" => 0)
# The paramets in order of the following input are n_events, Emin, Emax,
# volume (Dict, found above), and interaction_type ("cc", "nc", or "ccnc")
data, att = Event_Gen.generate_eventlist_cylinder(10, 1e19, 1e19, vol, "ccnc")
