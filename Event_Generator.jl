#!/usr/bin/julia

using Distributions
using Interpolations
using ArgParse

include("./Support_Funcs.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--N"
            help = "The number of neutrino events to be generated"
            arg_type = Int
            default = 10
        "--Emax"
            help = "The maximum neutrino energy (eV)."
            arg_type = Float64
            default = 1e19
        "--Emin"
            help = "The minimum neutrino energy (eV)"
            arg_type = Float64
            default = 1e17
        "--spectrum_type"
            help = "Defines probability distribution for which neutrino energies
            are generated"
            arg_type = String
            default = "log_uniform"
        "--volume"
            help = "Dictionary specifying simulation volume. Can be either a
            cylinder specified via keys
                *fiducial_rmin: float
                    lower r coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_rmax: float
                    upper r coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_zmin: float
                    lower z coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_zmax: float
                    upper z coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * full_rmin: float (optional)
                    lower r coordinate of simulated volume (if not set it is set to 1/3 of the fiducial volume, if second vertices are not activated it is set to the fiducial volume)
                * full_rmax: float (optional)
                    upper r coordinate of simulated volume (if not set it is set to the fiducial volume + the 95% quantile of the tau decay length, if second vertices are not activated it is set to the fiducial volume)
                * full_zmin: float (optional)
                    lower z coordinate of simulated volume (if not set it is set to the fiducial volume - the tau decay length, if second vertices are not activated it is set to the fiducial volume)
                * full_zmax: float (optional)
                    upper z coordinate of simulated volume (if not set it is set to 1/3 of the fiducial volume , if second vertices are not activated it is set to the fiducial volume)
            or a cube specified with
            * fiducial_xmin: float
                    lower x coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_xmax: float
                    upper x coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_ymin: float
                    lower y coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_ymax: float
                    upper y coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_zmin: float
                    lower z coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * fiducial_zmax: float
                    upper z coordinate of fiducial volume (the fiducial volume needs to be chosen large enough such that no events outside of it will trigger)
                * full_xmin: float (optional)
                    lower x coordinate of simulated volume (if not set it is set to the fiducial volume - the 95% quantile of the tau decay length, if second vertices are not activated it is set to the fiducial volume)
                * full_xmax: float (optional)
                    upper x coordinate of simulated volume (if not set it is set to the fiducial volume + the 95% quantile of the tau decay length, if second vertices are not activated it is set to the fiducial volume)
                * full_ymin: float (optional)
                    lower y coordinate of simulated volume (if not set it is set to the fiducial volume - the 95% quantile of the tau decay length, if second vertices are not activated it is set to the fiducial volume)
                * full_ymax: float (optional)
                    upper y coordinate of simulated volume (if not set it is set to the fiducial volume + the 95% quantile of the tau decay length, if second vertices are not activated it is set to the fiducial volume)
                * full_zmin: float (optional)
                    lower z coordinate of simulated volume (if not set it is set to 1/3 of the fiducial volume, if second vertices are not activated it is set to the fiducial volume)
                * full_zmax: float (optional)
                    upper z coordinate of simulated volume (if not set it is set to the fiducial volume - the tau decay length, if second vertices are not activated it is set to the fiducial volume)
                    "
            arg_type = Dict
            default = Dict("fiducial_rmin" => 0, "fiducial_rmax" => 5, "fiducial_zmin" => -2.7, "fiducial_zmax" => 0)
        "--output_name"
            help = "The name of the output file"
            arg_type = String
            default = "output.csv"
        "--thetamin"
            help = "Lower zenith angle for neutrino arrival direction (rad)"
            arg_type = Float64
            default = 0
        "--thetamax"
            help = "Upper zenith angle for neutrino arrival direction (rad)"
            arg_type = Float64
            default = 2*pi
        "--phimin"
            help = "Lower azimuth angle for neutrino arrival direction (rad)"
            arg_type = Float64
            default = 0
        "--phimax"
            help = "Upper azimuth angle for neutrino arrival direction (rad)"
            arg_type = Float64
            default = pi
        "--start_event-id"
            help = "Event number of first event"
            arg_type = Int64
            default = 1
        "--flavor"
            help = "Array of ints. Specify which neutrino flavors to generate.
            A uniform distribution of all specified flavors is assumed.
            *12: electron neutrino
            *14: muon neutrino
            *16: tau neutrino"
            arg_type = Array
            default = [12,-12,14,-14,16,-16]
        "--max_n_events_batch"
            help = "The maximum number of events that get generated per batch.
            Relevant if fiducial volume cut is applied."
            arg_type = Int64
            default = BigInt(1e6)
        "--write_events"
            help = "Choose whether to write results to a file."
            arg_type = Bool
            default = true
        "--interaction_type"
            help = "Interaction type. Default is 'ccnc' which randomly chooses
            neural current (NC) or charged-current (CC) interactions. User
            can also specify 'nc' or 'cc' to exclusively simulate NC or CC.
                * 'cc': charged-current
                * 'nc': neutral-current"
            arg_type = String
            default = "ccnc"
        "--start_event_id"
            help = "Event number of first event"
            arg_type = Int64
            default = 1
        end
    return parse_args(s)
end

"""
function generate_eventlist_cylinder(n_events, Emin, Emax, volume,
    thetamin=0, thetamax = pi, phimin=0, phimax=2*pi, start_event_id=1,
    flavor=[12,-12,14,-14,16,-16], spectrum="log_uniform",
    max_n_events_batch=1e5, write_events= true, interaction_type="cc")
"""

function generate_eventlist_cylinder()
    """
    Generates neutrino interactions (vertex positions, neutrino directions,
    neutrino flavor, charged/neutral current).
    """

    parsed_args = parse_commandLine()

    n_events = parsed_args["n_events"]
    Emin = parsed_args["Emin"]
    Emax = parsed_args["Emax"]
    volume = parsed_args["volume"]
    thetamin = parsed_args["thetamin"]
    thetamax = parsed_args["thetamax"]
    phimin = parsed_args["phimin"]
    phimax = parsed_args["phimax"]
    start_event_id = parsed_args["start_event_id"]
    flavor = parsed_args["flavor"]
    spectrum = parsed_args["spectrum"]
    max_n_events_batch = parsed_args["max_n_events_batch"]
    write_events = parsed_args["write_events"]
    interaction_type = parsed_args["interaction_type"]

    t_start = time()
    #attributes = Dict{String, Float64}()
    attributes = Dict()
    n_events = BigInt(n_events)

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
    n_batches = round(BigInt, ceil(n_events / max_n_events_batch))
    for i_batch in 1:n_batches
        data_sets = Dict()
        n_events_batch = round(BigInt, max_n_events_batch)

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
        data_sets["n_interaction"] = ones(BigInt, n_events_batch)
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
    print(attributes)
    return data_sets_fiducial, attributes
end

#vol = Dict("fiducial_rmin" => 0, "fiducial_rmax" => 5, "fiducial_zmin" => -2.7, "fiducial_zmax" => 0)
#data, att = generate_eventlist_cylinder(10, 1e18, 1e19, vol)

#Check what happens when Emax not greater than Emin
