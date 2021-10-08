
using Distributions

function set_volume_attributes(volume, attributes)
    """
    Interprets volume input.
    volume: dictionary
    proposal: bool
    attributes: dictionary
    """
    n_events = attributes['n_events']
    if (haskey(volume, 'fiducial_rmax')) #user specifies a cylinder
        if (haskey(volume, 'fiducial_rmin')
            attributes['fiducial_rmin'] = volume['fiducial_rmin']
        else
            attributes['fiducial_rmin'] = 0
        end

        attributes['fiducial_rmax'] = volume['fiducial_rmax']
        attributes['fiducial_zmin'] = volume['fiducial_zmin']
        attributes['fiducial_zmax'] = volume['fiducial_zmax']

        rmin = attributes['fiducial_rmin']
        rmax = attributes['fiducial_rmax']
        zmin = attributes['fiducial_zmin']
        zmax = attributes['fiducial_zmax']
        volume_fiducial = pi*(rmax^2 - rmin^2)*(zmax - zmin)

        if (haskey(volume, 'full_rmax'))
            rmax = volume['full_rmax']
        end
        if (haskey(volume, 'full_rmin'))
            rmin = volume['full_rmin']
        end
        if (haskey(volume, 'full_zmax'))
            zmax = volume['full_zmax']
        end
        if (haskey(volume, 'full_zmin'))
            zmin = volume['full_zmin']
        end

        volume_full = pi*(rmax^2 - rmin^2)*(zmax - zmin)
        # increase total # of events so we have same amount in fiducial vol
        n_events = round(Int32, n_events*volume_full / volume_fiducial)
        attributes['n_events'] = n_events

        attributes['rmin'] = rmin
        attributes['rmax'] = rmax
        attributes['zmin'] = zmin
        attributes['zmax'] = zmax

        V = pi*(rmax^2 - rmin^2)*(zmax - zmin)
        attributes['volume']  = V #save full sim vol to simplify eff vol calc
        attributes['area'] = pi*(rmax^2 - rmin^2)
    elseif (haskey(volume, 'fiducial_xmax')) #user specifies a cube
        attributes['fiducial_xmax'] = volume['fiducial_xmax']
        attributes['fiducial_xmin'] = volume['fiducial_xmin']
        attributes['fiducial_ymax'] = volume['fiducial_ymax']
        attributes['fiducial_ymin'] = volume['fiducial_ymin']
        attributes['fiducial_zmin'] = volume['fiducial_zmin']
        attributes['fiducial_zmax'] = volume['fiducial_zmax']

        xmin = attributes['fiducial_xmin']
        xmax = attributes['fiducial_xmax']
        ymin = attributes['fiducial_ymin']
        ymax = attributes['fiducial_ymax']
        zmin = attributes['fiducial_zmin']
        zmax = attributes['fiducial_zmax']
        volume_fiducial = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        if (haskey(volume, 'full_xmax'))
            xmin = volume['full_xmin']
            xmax = volume['full_xmax']
            ymin = volume['full_ymin']
            ymax = volume['full_ymax']
            zmin = volume['full_zmin']
            zmax = volume['full_zmax']
        end

        volume_full = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        n_events = round(Int32, n_events * volume_full / volume_fiducial)
        attributes['n_events'] = n_events

        attributes['xmin'] = xmin
        attributes['xmax'] = xmax
        attributes['ymin'] = ymin
        attributes['ymax'] = ymax
        attributes['zmin'] = zmin
        attributes['zmax'] = zmax

        V = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        attributes['volume'] = V  # save full sim vol to simplify eff vol calc
        attributes['area'] = (xmax - xmin) * (ymax - ymin)
    else
        throw(DomainError())
    end
end

function generate_vertex_positions(attributes, n_events, rnd=nothing)
    """
    Generates vertex positions randomly distributed in simulation volume
    and outputs relevent quantities.
    """
    if (haskey(attributes, 'fiducial_rmax'))
        rr_full = rand(Uniform())


end

function generate_eventlist_cylinder(filename, n_events, Emin, Emax, volume,
    thetamin=0, thetamax = pi, phimin=0, phimax=2*pi, start_event_id=1,
    flavor=[12,-12,14,-14,16,-16], n_events_per_file=nothing,
    spectrum='log_uniform', start_file_id=0, max_n_events_batch=1e5,
    write_events= true, seed=nothing, interaction_type='ccnc')

    """
    Generates neutrino interactions (vertex positions, neutrino directions,
    neutrino flavor, charged/neutral current).
    """

    t_start = time()
    attributes = Dict{String, Float64}()

    attributes['start_event_id'] = start_event_id
    attributes['n_events'] = n_events
    attributes['flavors'] = flavor
    attributes['Emin'] = Emin
    attributes['Emax'] = Emax
    attributes['thetamin'] = thetamin
    attributes['thetamax'] = thetamax
    attributes['phimin'] = phimin
    attributes['phimax'] = phimax

    data_sets = Dict()
    data_sets_fiducial = Dict()

    time_proposal = 0

    set_volume_attributes(volume, attributes=attributes)

end
