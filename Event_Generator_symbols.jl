#!/usr/bin/julia
module Event_Gen_sym
export generate_eventlist_cylinder

using Revise
using Distributions
using Interpolations
using ArgParse
using CSV, Tables
using DataFrames
using HDF5
using Distributed
using FLoops

include("./Support_Funcs_symbols.jl")

"""
Generates neutrino interactions (vertex positions, neutrino directions,
neutrino flavor, charged/neutral current).
"""

function generate_eventlist_cylinder(n_events, Emin::Float64,
    Emax::Float64, volume, interaction_type, thetamin=0.0, thetamax = 1*pi,
    phimin=0, phimax=2*pi, start_event_id=1, flavor=[12,-12,14,-14,16,-16],
    n_events_per_file=nothing, start_file_id=0, max_n_events_batch=256,
    write_events= true)

    attributes = Dict()

    attributes[:start_event_id] = start_event_id::Int
    attributes[:n_events] = Int(n_events) #ensure type stability
    attributes[:flavors] = flavor::Vector{Int64}
    attributes[:Emin] = Emin::Float64
    attributes[:Emax] = Emax::Float64
    attributes[:thetamin] = thetamin
    attributes[:thetamax] = thetamax
    attributes[:phimin] = phimin
    attributes[:phimax] = phimax


    data_sets = Dict()
    data_sets_fiducial = Dict()

    #generate detector volume attributes
    attributes = set_volume_attributes(volume, attributes)
    n_events = attributes[:n_events]

    #split number of events into batches
    n_batches = round(Int, ceil(n_events / max_n_events_batch))

    #generate event attributes in batches
    @floop for i_batch in 0:(n_batches - 1)
        data_sets = Dict()
        n_events_batch = round(Int, max_n_events_batch)

        if (i_batch + 1 == n_batches) #check if this is the last batch
            n_events_batch = round(Int, n_events - (i_batch * max_n_events_batch))
        end

        #generate vertex positions in cartesian coordinates
        data_sets[:xx], data_sets[:yy], data_sets[:zz] = generate_vertex_positions(attributes, n_events_batch)
        #data_sets[:zz] = zero(data_sets[:zz]) #muons interact at the surface so zz => 0

        # generate neutrino vertices randomly
        data_sets[:azimuths] = rand(Uniform(phimin, phimax), n_events_batch)
        # zenith directions are distributed as sin(theta) (to make dist. isotropic) * cos(theta) (to acc for projection onto surf)
        data_sets[:zeniths] = acos.(rand(Uniform(cos(thetamax), cos(thetamin)), n_events_batch)).^0.5

        #label each event with an ID
        data_sets[:event_group_ids] = collect((i_batch*max_n_events_batch):((i_batch*max_n_events_batch)+n_events_batch)-1).+start_event_id
        #count number of interactions(?)
        data_sets[:n_interaction] = ones(Int, n_events_batch)
        data_sets[:vertex_times] = zeros(Float64, n_events_batch)

        #generate neutrino flavors randomly
        p = repeat([1/length(flavor)], length(flavor)) #create probability vector for flavors
        rng = Categorical(p)
        data_sets[:flavors] = flavor[rand(rng, n_events_batch)]

        #generate neutrino energies randomly
        data_sets[:energies] = get_energies(n_events_batch, Emin, Emax)

        #generate charged/neutral current randomly
        if interaction_type == "ccnc"
            data_sets[:interaction_type] = get_ccnc(n_events_batch)
        elseif interaction_type == "cc"
            data_sets[:interaction_type] = repeat([:cc], outer=[n_events_batch])
        elseif interaction_type == "nc"
            data_sets[:interaction_type] = repeat([:nc], outer=[n_events_batch])
        end

        #generate inelasticity (fraction of neutrino energy that goes into shower)
        data_sets[:inelasticity] = get_neutrino_inelasticity(n_events_batch)
        #get true energy by multiplying original E by inelasticity
        data_sets[:shower_energies] = data_sets[:energies] .* data_sets[:inelasticity]

        #set all neutrino interaction results to hadronic showers by default
        data_sets[:shower_type] = repeat(["had"], n_events_batch)

        # now add EM showers when appropriate
        #check which showers are electromagnetic by requiring an electron neutrino (12) and cc interaction
        em_shower_mask = (data_sets[:interaction_type] .== "cc").*(abs.(data_sets[:flavors]) .== 12)
        n_inserted = 0

        for i in keepat!(collect(0:(n_events_batch - 1)), em_shower_mask)  # loop over all events where an EM shower needs to be inserted
            for (key, value) in data_sets
                insert!(data_sets[key], (i+1) + 1 + n_inserted, data_sets[key][(i+1) + n_inserted])  # copy event
            end
            data_sets[:shower_energies][i + 1 + n_inserted] = (1 - data_sets[:inelasticity][i + 1 + n_inserted]) * data_sets[:energies][i + 1 + n_inserted]

			data_sets[:shower_type][i + 1 + n_inserted] = "em"
            n_inserted += 1
        end



        if n_batches == 1
            data_sets_fiducial = data_sets
        else
            for (key, value) in data_sets
                if key ∉ keys(data_sets_fiducial)
                    data_sets_fiducial[key] = data_sets[key]
                else
                    vcat(data_sets_fiducial[key], data_sets[key])
                end
            end
        end

    end

    # make event group ids consecutive - useful if secondary interactions simulated
    # where many of the initially generated neutrinos don't end up in fiducial vol
    for (key, value) in data_sets
        if key ∉ keys(data_sets_fiducial)
            data_sets_fiducial[key] = data_sets[key]
        end
    end
    egids = data_sets_fiducial[:event_group_ids]
    uegids = unique(egids)
    uegids_inverse = []
    for value in egids
        if value in uegids
            append!(uegids_inverse, findall(x -> x == value, uegids)[1])
        end
    end

    # assign every shower a unique ID
    data_sets_fiducial[:shower_ids] = collect(range(0, step=1, length(data_sets_fiducial[:shower_energies])))
    # make the event group ids consecutive, this is useful if secondary interactions are simulated where many of the
    # initially generated neutrinos don't end up in the fiducial volume
    data_sets_fiducial[:event_group_ids] = uegids_inverse .+ start_event_id

    """
    NOT CURRENTLY WORKING:
    if write_events
        write_events_to_hdf5(filename="testing", data_sets=data_sets_fiducial,
            attributes=attributes, n_events_per_file=n_events_per_file,
            start_file_id=start_file_id)
    end
    """

    CSV.write("data_output.csv", data_sets_fiducial, header=false)
    CSV.write("attributes_output.csv", attributes, header=false)
end

vol = Dict(:fiducial_rmin => 0.0, :fiducial_rmax => 5000.0, :fiducial_zmin => -2700, :fiducial_zmax => 0.0)

function benchmark_my_function(f::Function, num_list::Vector{T} ) where { T<:Integer }
	times_list = zeros(length(num_list))
	for (i,n) in enumerate(num_list)
		times_list[i] = @belapsed $f($n, 1e18, 1e18, vol, "ccnc")
	end
	return (;num_list, times_list)
end


end
