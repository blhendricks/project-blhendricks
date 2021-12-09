using Distributions
using Interpolations
using ArgParse

"""
Get nature of interaction current: cc or nc
"""
function get_ccnc(n_events)
    random_sequence = rand(Uniform(0.0, 1.0), n_events)
    ccnc = fill("nc", 1, n_events)

    for (i, r) in enumerate(random_sequence)
        if(r <= 0.7064)
            ccnc[i] = "cc"
        end
    end

    return vec(ccnc) #cast as vector to avoid breaking other code
end

"""
Standard inelasticity for deep inelastic scattering
"""
function get_neutrino_inelasticity(n_events)
    R1 = 0.36787944
    R2 = 0.63212056
    inelasticities = (-log.(R1 .+ (rand(Uniform(0.0, 1.0), n_events).*R2))).^2.5
    return inelasticities
end

"""
Generates a random distribution of energies following a certain spectrum

Params
Emin, Emax: float
n_event: int
flux: function
"""
function get_energies(n_events, Emin, Emax)
    if Emin == Emax
            energies = 10 .^ (repeat([log10(Emin)], n_events))
    else
    energies = 10 .^ (rand(Uniform(log10(Emin), log10(Emax)), n_events))
    end

    return energies
end

"""
Interprets volume input.
volume: dictionary
proposal: bool
attributes: dictionary
"""
function set_volume_attributes(volume, attributes)
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

"""
Generates vertex positions randomly distributed in simulation volume
and outputs relevent quantities.
"""
function generate_vertex_positions(attributes, n_events)
    if (haskey(attributes, "fiducial_rmax"))
        print("rmin squared = ", attributes["rmin"]^2, "\n")
        print("rmax squared = ", attributes["rmax"]^2, "\n")
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
