
using Distributions
using Interpolations

function get_neutrino_inelasticity(n_events)
    """
    Standard inelasticity for deep inelastic scattering
    """
    R1 = 0.36787944
    R2 = 0.63212056
    inelasticities = (-log.(R1 .+ (rand(Uniform(0.0, 1.0), n_events).*R2))).^2.5
    return inelasticities
end

function get_energy_from_flux(Emin, Emax, n_events, flux)
    xx_edges = collect(range(Emin, Emax, 10000000))
    xx = 0.5 * (xx_edges[2:end] .+ xx_edges[1:(end-1)])
    yy = flux(xx)
    cum_values = zeros(size(xx_edges))
    cum_values[2:end] = cumsum(yy * diff(xx_edges))
    inv_cdf = LinearInterpolation(cum_values, xx_edges)
    r = rand(Uniform(0, maximum(cum_values)), n_events)
    return inv_cdf(r)
end

function get_energies(n_events, Emin, Emax, spectrum_type)
    """
    Generates a random distribution of energies following a certain spectrum

    Params
    Emin, Emax: float
    n_event: int
    flux: function
    """

    if spectrum_type == "log_uniform"
        energies = 10 .^ (rand(Uniform(log10(Emin), log10(Emax)), n_events))
    elseif startswith(spectrum_type, "E-") # generate an E^gamma spectrum
        gamma = float(spectrum_type[2:end])
        gamma += 1
        Nmin = (Emin)^gamma
        Nmax = (Emax)^gamma

        function get_inverse_spectrum(N, gamma)
            return exp(log(N)/gamma)
        end

        energies = get_inverse_spectrum(rand(Uniform(Nmax, Nmin), n_events), gamma)

    elseif spectrum_type == "GZK-1"
        energies = get_energy_from_flux(Emin, Emax, n_events, get_GZK_1)
    elseif spectrum_type == "IceCube-nu-2017"
        energies = get_energy_from_flux(Emin, Emax, n_events, ice_cube_nu_fit)
    elseif spectrum_type == "GZK-1+IceCube-nu-2017"

        function J(E)
            return ice_cube_nu_fit(E) + get_GZK_1(E)
        end

        energies = get_energy_from_flux(Emin, Emax, n_events, J)
    else
        println("Passed spectrum not implemented.")
        throw(DomainError())
        #throw("unimplemented") more appropriate?
    end
    return energies
end

function set_volume_attributes(volume, attributes)
    """
    Interprets volume input.
    volume: dictionary
    proposal: bool
    attributes: dictionary
    """
    n_events = attributes["n_events"]

    if (haskey(volume, "fiducial_rmax")) #user specifies a cylinder
        if (haskey(volume, "fiducial_rmin"))
            attributes["fiducial_rmin"] = volume["fiducial_rmin"]
        else
            attributes["fiducial_rmin"] = 0
        end

        attributes["fiducial_rmax"] = volume["fiducial_rmax"]
        attributes["fiducial_zmin"] = volume["fiducial_zmin"]
        attributes["fiducial_zmax"] = volume["fiducial_zmax"]

        rmin = attributes["fiducial_rmin"]
        rmax = attributes["fiducial_rmax"]
        zmin = attributes["fiducial_zmin"]
        zmax = attributes["fiducial_zmax"]
        volume_fiducial = pi*(rmax^2 - rmin^2)*(zmax - zmin)

        if (haskey(volume, "full_rmax"))
            rmax = volume["full_rmax"]
        end
        if (haskey(volume, "full_rmin"))
            rmin = volume["full_rmin"]
        end
        if (haskey(volume, "full_zmax"))
            zmax = volume["full_zmax"]
        end
        if (haskey(volume, "full_zmin"))
            zmin = volume["full_zmin"]
        end

        volume_full = pi*(rmax^2 - rmin^2)*(zmax - zmin)
        # increase total # of events so we have same amount in fiducial vol
        n_events = round(BigInt, n_events*volume_full / volume_fiducial)
        attributes["n_events"] = n_events

        attributes["rmin"] = rmin
        attributes["rmax"] = rmax
        attributes["zmin"] = zmin
        attributes["zmax"] = zmax

        V = pi*(rmax^2 - rmin^2)*(zmax - zmin)
        attributes["volume"]  = V #save full sim vol to simplify eff vol calc
        attributes["area"] = pi*(rmax^2 - rmin^2)
    elseif (haskey(volume, "fiducial_xmax")) #user specifies a cube
        attributes["fiducial_xmax"] = volume["fiducial_xmax"]
        attributes["fiducial_xmin"] = volume["fiducial_xmin"]
        attributes["fiducial_ymax"] = volume["fiducial_ymax"]
        attributes["fiducial_ymin"] = volume["fiducial_ymin"]
        attributes["fiducial_zmin"] = volume["fiducial_zmin"]
        attributes["fiducial_zmax"] = volume["fiducial_zmax"]

        xmin = attributes["fiducial_xmin"]
        xmax = attributes["fiducial_xmax"]
        ymin = attributes["fiducial_ymin"]
        ymax = attributes["fiducial_ymax"]
        zmin = attributes["fiducial_zmin"]
        zmax = attributes["fiducial_zmax"]
        volume_fiducial = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        if (haskey(volume, "full_xmax"))
            xmin = volume["full_xmin"]
            xmax = volume["full_xmax"]
            ymin = volume["full_ymin"]
            ymax = volume["full_ymax"]
            zmin = volume["full_zmin"]
            zmax = volume["full_zmax"]
        end

        volume_full = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        n_events = round(BigInt, n_events * volume_full / volume_fiducial)
        attributes["n_events"] = n_events

        attributes["xmin"] = xmin
        attributes["xmax"] = xmax
        attributes["ymin"] = ymin
        attributes["ymax"] = ymax
        attributes["zmin"] = zmin
        attributes["zmax"] = zmax

        V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        attributes["volume"] = V  # save full sim vol to simplify eff vol calc
        attributes["area"] = (xmax - xmin) * (ymax - ymin)
    else
        println("'fiducial_rmin' or 'fiducial_rmax' not part of 'attributes'")
        throw(DomainError())
    end

    return attributes
end

function generate_vertex_positions(attributes, n_events)
    """
    Generates vertex positions randomly distributed in simulation volume
    and outputs relevent quantities.
    """

    if (haskey(attributes, "fiducial_rmax"))
        rr_full = rand(Uniform(attributes["rmin"]^2, attributes["rmax"]^2), n_events).^0.5
        phiphi = rand(Uniform(0, 2*pi), n_events)
        xx = rr_full .* cos.(phiphi)
        yy = rr_full .* sin.(phiphi)
        zz = rand(Uniform(attributes["zmin"], attributes["zmax"]), n_events)
        return xx, yy, zz
    elseif (haskey(attributes, "fiducial_xmax"))
        xx = rand(Uniform(attributes["xmin"], attributes["xmax"]), n_events)
        yy = rand(Uniform(attributes["ymin"], attributes["ymax"]), n_events)
        zz = rand(Uniform(attributes["zmin"], attributes["zmax"]), n_events)
        return xx, yy, zz
    else
        println("'fiducial_rmin' or 'fiducial_rmax' not part of 'attributes'")
        throw(DomainError())
    end
end

#generate_vertex_positions(att, n_events)

function generate_eventlist_cylinder(n_events, Emin, Emax, volume,
    thetamin=0, thetamax = pi, phimin=0, phimax=2*pi, start_event_id=1,
    flavor=[12,-12,14,-14,16,-16], n_events_per_file=nothing,
    spectrum="log_uniform", start_file_id=0, max_n_events_batch=1e5,
    write_events= true, seed=nothing, interaction_type="ccnc")

    """
    Generates neutrino interactions (vertex positions, neutrino directions,
    neutrino flavor, charged/neutral current).
    """

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
        """ Inelasticities library not converted to Julia yet
        if interaction_type == "ccnc"
            data_sets["interaction_type"] = inelasticities.get_ccnc(n_events_batch)
        """
        if interaction_type == "cc"
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

    return data_sets_fiducial, attributes
end

vol = Dict("fiducial_rmin" => 0, "fiducial_rmax" => 5, "fiducial_zmin" => -2.7, "fiducial_zmax" => 0)
data, att = generate_eventlist_cylinder(1e10, 1e18, 1e19, vol)

#Check what happens when Emax not greater than Emin
