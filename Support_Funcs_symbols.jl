using Distributions
using Interpolations
using ArgParse
using Printf
using HDF5
using BenchmarkTools

"""
Get nature of interaction current: cc or nc
"""
function get_ccnc(n_events::Int)
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
function get_neutrino_inelasticity(n_events::Int)
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
function get_energies(n_events::Int, Emin::Float64, Emax::Float64)
    if Emin == Emax
            energies = 10.0 .^ (repeat([log10(Emin)], n_events))
    else
    energies = 10.0 .^ (rand(Uniform(log10(Emin), log10(Emax)), n_events))
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
    n_events = attributes[:n_events]
    if (haskey(volume, :fiducial_rmax)) #user specifies a cylinder
        if (haskey(volume, :fiducial_rmin))
            attributes[:fiducial_rmin] = volume[:fiducial_rmin]
        else
            attributes[:fiducial_rmin] = 0.0
        end

        attributes[:fiducial_rmax] = volume[:fiducial_rmax]
        attributes[:fiducial_zmin] = volume[:fiducial_zmin]
        attributes[:fiducial_zmax] = volume[:fiducial_zmax]

        rmin = float(attributes[:fiducial_rmin])
        rmax = float(attributes[:fiducial_rmax])
        zmin = float(attributes[:fiducial_zmin])
        zmax = float(attributes[:fiducial_zmax])
        volume_fiducial = pi*(rmax^2 - rmin^2)*(zmax - zmin)

        if (haskey(volume, :full_rmax))
            rmax = float(volume[:full_rmax])
        end
        if (haskey(volume, :full_rmin))
            rmin = float(volume[:full_rmin])
        end
        if (haskey(volume, :full_zmax))
            zmax = float(volume[:full_zmax])
        end
        if (haskey(volume, :full_zmin))
            zmin = float(volume[:full_zmin])
        end

        volume_full = pi*(rmax^2 - rmin^2)*(zmax - zmin)
        # increase total # of events so we have same amount in fiducial vol
        n_events = round(BigInt, n_events*volume_full / volume_fiducial)
        attributes[:n_events] = n_events

        attributes[:rmin] = rmin
        attributes[:rmax] = rmax
        attributes[:zmin] = zmin
        attributes[:zmax] = zmax

        V = pi*(rmax^2 - rmin^2)*(zmax - zmin)
        attributes[:volume]  = V #save full sim vol to simplify eff vol calc
        attributes[:area] = pi*(rmax^2 - rmin^2)
    elseif (haskey(volume, :fiducial_xmax)) #user specifies a cube
        attributes[:fiducial_xmax] = float(volume[:fiducial_xmax])
        attributes[:fiducial_xmin] = float(volume[:fiducial_xmin])
        attributes[:fiducial_ymax] = float(volume[:fiducial_ymax])
        attributes[:fiducial_ymin] = float(volume[:fiducial_ymin])
        attributes[:fiducial_zmin] = float(volume[:fiducial_zmin])
        attributes[:fiducial_zmax] = float(volume[:fiducial_zmax])

        xmin = attributes[:fiducial_xmin]
        xmax = attributes[:fiducial_xmax]
        ymin = attributes[:fiducial_ymin]
        ymax = attributes[:fiducial_ymax]
        zmin = attributes[:fiducial_zmin]
        zmax = attributes[:fiducial_zmax]
        volume_fiducial = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        if (haskey(volume, :full_xmax))
            xmin = float(volume[:full_xmin])
            xmax = float(volume[:full_xmax])
            ymin = float(volume[:full_ymin])
            ymax = float(volume[:full_ymax])
            zmin = float(volume[:full_zmin])
            zmax = float(volume[:full_zmax])
        end

        volume_full = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        n_events = round(BigInt, n_events * volume_full / volume_fiducial)
        attributes[:n_events] = n_events

        attributes[:xmin] = xmin
        attributes[:xmax] = xmax
        attributes[:ymin] = ymin
        attributes[:ymax] = ymax
        attributes[:zmin] = zmin
        attributes[:zmax] = zmax

        V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        attributes[:volume] = V  # save full sim vol to simplify eff vol calc
        attributes[:area] = (xmax - xmin) * (ymax - ymin)
    else
        println("'fiducial_rmin' or 'fiducial_rmax' not part of 'attributes'")
    end

    return attributes
end

"""
Generates vertex positions randomly distributed in simulation volume
and outputs relevent quantities.
"""
function generate_vertex_positions(attributes, n_events::Int)
    if (haskey(attributes, :fiducial_rmax))
        rr_full = rand(Uniform(attributes[:rmin]^2, attributes[:rmax]^2), n_events).^0.5
        phiphi = rand(Uniform(0, 2*pi), n_events)
        xx = rr_full .* cos.(phiphi)
        yy = rr_full .* sin.(phiphi)
        zz = rand(Uniform(attributes[:zmin], attributes[:zmax]), n_events)
        return xx, yy, zz
    elseif (haskey(attributes, :fiducial_xmax))
        xx = rand(Uniform(attributes[:xmin], attributes[:xmax]), n_events)
        yy = rand(Uniform(attributes[:ymin], attributes[:ymax]), n_events)
        zz = rand(Uniform(attributes[:zmin], attributes[:zmax]), n_events)
        return xx, yy, zz
    else
        println("'fiducial_rmin' or 'fiducial_rmax' not part of 'attributes'")
    end
end


"""
Writes NuRadioMC data to hdf5 file
Can automatically split dataset up into multiple files for multiprocessing
"""

"""
function write_events_to_hdf5(filename, data_sets, attributes, n_events_per_file=nothing,
    start_file_id=0)

    print("test\n\n\n")

    n_events = attributes[:n_events]
    total_number_of_events = attributes[:n_events]

    if "start_event_id" ∉ keys(attributes)
        attributes[:start_event_id] = 0
    end

    if n_events_per_file == nothing
        n_events_per_file = n_events
    else
        n_events_per_file = round(Int, n_events_per_file)
    end
    iFile = -1
    evt_id_first = data_sets[:event_group_ids][1]
    evt_id_last_previous = 0 #save last ID of previous file
    start_index = 0
    n_events_total = 0

    while true
        iFile += 1
        filename2 = filename
        evt_ids_this_file = unique(data_sets[:event_group_ids])[iFile * n_events_per_file: (iFile+1)*n_events_per_file]
        if(length(evt_ids_this_file) == 0)
            break #no more events to write in file
        end

        if((iFile > 0) || (n_events_per_file < n_events))
            filename2 = filename*@sprintf(".part%04d", iFile + start_file_id)
        end
        fout = h5open(filename2, "w")

        attributes(fout)[:VERSION MAJOR] = 2
        attributes(fout)[:VERSION MAJOR] = 2
        attributes(fout)[:header] = """
# all quantities are in the default NuRadioMC units (i.e., meters, radians and eV)
# all geometry quantities are in the NuRadioMC default local coordinate system:
#     coordinate origin is at the surface
#     x axis is towards Easting, y axis towards Northing, z axis upwards
#     zenith/theta angle is defined with respect to z axis, i.e. 0deg = upwards, 90deg = towards horizon, 180deg = downwards
#     azimuth/phi angle counting northwards from East
#
# the collumns are defined as follows
# 1. event id (integer)
# 2. neutrino flavor (integer) encoded as using PDG numbering scheme, particles have positive sign, anti-particles have negative sign, relevant for us are:
#       12: electron neutrino
#       14: muon neutrino
#       16: tau neutrino
# 3. energy of neutrino (double)
# 4. charge or neutral current interaction (string, one of ['cc', 'nc']
# 5./6./7. position of neutrino interaction vertex in cartesian coordinates (x, y, z) (in default NuRadioMC local coordinate system)
# 8. zenith/theta angle of neutrino direction (pointing to where it came from, i.e. opposite to the direction of propagation)
# 9. azimuth/phi angle of neutrino direction (pointing to where it came from, i.e. opposite to the direction of propagation)
# 10. inelasticity (the fraction of neutrino energy that goes into the hadronic part)
#
"""
        for (key, value) in collect(attributes)
            attributes(fout)[key] = value
        end
        attributes(fout)[:total_number_of_events] = total_number_of_events

        evt_id_first = evt_ids_this_file[1]
        evt_id_last = last(evt_ids_this_file)

        tmp = dropdims(findall(x->x==evt_id_last, data_sets[:event_group_ids]), tuple(findall(size(a) .== 1)...))
        if size(tmp) == (1,)
            stop_index = last(tmp) + 1
        end

        for (key, value) in collect(data_sets)
            fout[key] = value[start_index: stop_index]
        end

        evt_ids_next_file =unique(data_sets[:event_group_ids])[(iFile+1) * n_events_per_file: (iFile+2) * n_events_per_file]
        n_events_this_file = nothing
        if (iFile == 0 && length(evt_ids_next_file) == 0)
            n_events_this_file = total_number_of_events
        elseif (length(evt_ids_next_file) == 0)
            n_events_this_file = total_number_of_events - (evt_id_last_previous + 1) + attributes[:start_event_id]
        elseif (iFile == 0)
            n_events_this_file = evt_id_last - attributes[:start_event_id] + 1
        else
            n_events_this_file = evt_id_last - evt_id_last_previous
        end

        attributes(fout)[:n_events] = n_events_this_file
        fout.close()
        n_events_total += n_events_this_file

        start_index = stop_index

        evt_id_last_previous = evt_id_last
        if (evt_id_last == n_events)
            break
        end
    end
end
"""
