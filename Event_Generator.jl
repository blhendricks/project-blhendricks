#!/usr/bin/julia

using Distributions
using Interpolations
using ArgParse
using CSV, Tables
using DataFrames

include("./Support_Funcs.jl")


function generate_eventlist_cylinder(n_events, Emin, Emax, volume,
    interaction_type, thetamin=0, thetamax = pi, phimin=0, phimax=2*pi,
    start_event_id=1, flavor=[12,-12,14,-14,16,-16], spectrum="log_uniform",
    max_n_events_batch=1e5, write_events= true)

    """
    Generates neutrino interactions (vertex positions, neutrino directions,
    neutrino flavor, charged/neutral current).
    """

    t_start = time()
    #attributes = Dict{String, Float64}()
    attributes = Dict()
    n_events = Int(n_events)

    attributes["start_event_id"] = start_event_id
    attributes["n_events"] = n_events
    attributes["flavors"] = flavor
    attributes["Emin"] = Emin
    attributes["Emax"] = Emax
    attributes["thetamin"] = thetamin
    attributes["thetamax"] = thetamax
    attributes["phimin"] = phimin
    attributes["phimax"] = phimax

    data_sets = Dict()
    data_sets_fiducial = Dict()

    time_proposal = 0

    #n_events_batch = round(Int32, n_events_batch)

    attributes = set_volume_attributes(volume, attributes)
    n_events = attributes["n_events"]
    n_batches = round(Int, ceil(n_events / max_n_events_batch))
    for i_batch in 1:n_batches
        data_sets = Dict()
        n_events_batch = round(Int, max_n_events_batch)

        if (i_batch + 1) == n_batches
            n_events_batch = n_events - (i_batch * max_n_events_batch)
        end
        data_sets["xx"], data_sets["yy"], data_sets["zz"] = generate_vertex_positions(attributes, n_events_batch)
        data_sets["zz"] = zero(data_sets["zz"]) #muons interact at the surface so zz => 0

        # generate neutrino vertices randomly
        data_sets["azimuths"] = rand(Uniform(phimin, phimax), n_events_batch)
        # zenith directions are distributed as sin(theta) (to make dist. isotropic) * cos(theta) (to acc for projection onto surf)
        data_sets["zeniths"] = asin.(rand(Uniform(sin(thetamin)^2, sin(thetamax)^2), n_events_batch)).^0.5

        data_sets["event_group_ids"] = collect((i_batch*max_n_events_batch):((i_batch*max_n_events_batch)+n_events_batch-1)).+start_event_id
        data_sets["n_interaction"] = ones(Int, n_events_batch)
        data_sets["vertex_times"] = zeros(Float64, n_events_batch)

        #generate neutrino flavors randomly
        data_sets["flavors"] = flavor[rand(1:end, n_events_batch)]
        #generate neutrino energies randomly
        data_sets["energies"] = get_energies(n_events_batch, Emin, Emax, spectrum)
        #generate charged/neutral current randomly
        if interaction_type == "ccnc"
            data_sets["interaction_type"] = get_ccnc(n_events_batch)
        elseif interaction_type == "cc"
            data_sets["interaction_type"] = repeat(["cc"], outer=[n_events_batch])
        elseif interaction_type == "nc"
            data_sets["interaction_type"] = repeat(["nc"], outer=[n_events_batch])
        end

        #generate inelasticity
        data_sets["inelasticity"] = get_neutrino_inelasticity(n_events_batch)
        """
        This EM shower portion will be added later

        #add EM showers if appropriate
        em_shower_mask = (data_sets["interaction_type"] == "cc") & (abs.(data_sets["flavors"]) == 12)

        n_inserted = 0
        if em_shower_mask == true #loop over all events where EM shower needs to be inserted
            for i in collect(range(1, step=1, n_events_batch))
                for key in data_sets
                    data_sets[key]
        """
    end

    time_per_evt = time_proposal / (n_events + 1)

    # assign every shower a unique ID
    """
    data_sets_fiducial["shower_ids"] = collect(range(0, step=1, length(data_sets_fiducial["shower_energies"])))
    """
    # make event group ids consecutive - useful if secondary interactions simulated
    # where many of the initially generated neutrinos don't end up in fiducial vol
    for (key, value) in data_sets
        if key ∉ keys(data_sets_fiducial)
            data_sets_fiducial[key] = data_sets[key]
        end
    end
    egids = data_sets_fiducial["event_group_ids"]
    uegids = unique(egids)
    uegids_inverse = []
    for value in egids
        if value in uegids
            append!(uegids_inverse, findall(x -> x == value, uegids)[1])
        end
    end

    data_sets_fiducial["event_group_ids"] = uegids_inverse .+ start_event_id
    """
    Will put in option to write to file later
    if(write_events): ...

    for (key, value) in
    """

    CSV.write("data_output.csv", data_sets_fiducial, header=false)
    CSV.write("attributes_output.csv", attributes, header=false)
end

vol = Dict("fiducial_rmin" => 0, "fiducial_rmax" => 5, "fiducial_zmin" => -2.7, "fiducial_zmax" => 0)
data, att = generate_eventlist_cylinder(n_events=10, Emin=1e18, Emax=1e19,
    volume=vol, interaction_type="ccnc")
