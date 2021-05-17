module AtomDist

using Distributions

function initial_distribution(z_velocity::Float64,z_spread::Float64,
                              x_velocity::Float64,x_spread::Float64,
                              excitation_time::Float64,pi_time::Float64,
                              angle_offset::Float64,number_atoms::Int64,
                              laser_x_position::Float64
                              )
    """
    Returns a random sample of atoms meeting the distribution.

    Parameters
    =========
    z_velocity::Float64
        the average z velocity
    z_spread::Float64
        the std on z velocity
    x_velocity::Float64
        the average x velocity
    x_spread::Float64
            the std on x velocity
    excitation_time::Float64
        the duration for which atoms are being excited
    pi_time::Float64
        the time between start of exciting till sequence start.
    angle_offset::Float64
        how misaligned the plates are with the beam.
    number_atoms::Int64
        the number of atoms to produce.
    laser_x_position::Float64
        how far along the laser is

    Returns
    ========
    velz_x : Array{Float64,1d}
        distribution of velocities in x direction
    velz_z : Array{Float64,1d}
        distribution of velocities in z direction
    posz_z : Array{Float64,1d}
        distribution of positions in x direction
    posz_z : Array{Float64,1d}
        distribution of positions in x direction
    """

    velz = get_initial_velocities(number_atoms,x_velocity,x_spread,z_velocity,z_spread)
    velz_x = velz[1]
    velz_z = velz[2]
    velz = rotate_into_frame(velz_x,velz_z,angle_offset)
    velz_x = velz[1]
    velz_z = velz[2]
    posz = get_initial_positions(number_atoms,velz_x,velz_z,excitation_time,pi_time,laser_x_position)
    posz_x=  posz[1]
    posz_z = posz[2]
    return (velz_x,velz_z,posz_x,posz_z)
end


function get_initial_velocities(number_atoms::Int64,x_velocity::Float64,x_spread::Float64,z_velocity::Float64,z_spread::Float64)
    """
    Samples a gaussian to get a velocity distribution in the lab frame. returns different distributions in the x and z directions

    Parameters
    ===========
    number_atoms
        total number of atoms in simulation
    x_velocity
        mean speed in x direction
    x_spread
        standard deviation in x direction
    z_velocity
        the mean speed in z direction
    z_spread
        standard deviation in z direction

    Return
    =======
    (x_vels,z_vels) (Array{Float64,1d},Array{Float64,1d})
        returns a tuple containing two 1d vectors of the velocities
    """
    x_vels = rand(Normal(x_velocity,x_spread),number_atoms)
    z_vels = rand(Normal(z_velocity,z_spread),number_atoms)

    return (x_vels,z_vels)
end

function get_initial_positions(number_atoms::Int64,
                        initial_x_velocities::Array{Float64,1},initial_z_velocities::Array{Float64,1},
                        excitation_time::Float64,pi_time::Float64,laser_x_position::Float64)
    """
    Gets the x and z positions of the initial distribution. samples a uniform distribution between
    excitation times up to excitation time. computes the amount of time atoms have left of travel
    before the microwave pulse is applied and multiplies this by the atoms initial velocity.
    """
    excitation_array = rand(Uniform(0.0, excitation_time),number_atoms)
    journey_time = pi_time .- excitation_array
    x_position_array = journey_time.*initial_x_velocities .+laser_x_position
    z_position_array =  journey_time.*initial_z_velocities

    return (x_position_array,z_position_array)
end

function rotate_into_frame(initial_x_velocity::Array{Float64,1},initial_z_velocity::Array{Float64,1},angle_offset::Float64)
    """
    transforms the velocities into a rotated frame about the y axis. Into order to account for misalignments
    in the electrodes relative to the supersonic beam.

    Parameters
    ==========
    initial_x_velocity
        array of all atoms initial velocities in x
    initial_z_velocity
        array of all atoms initial velocities in z
    angle_offset
        the angle by which to rotate (radians)

    Returns
    (x_vels,z_vels) (Array{Float64,1d},Array{Float64,1d})
        returns a tuple containing the new 1d vectors of the velocities
    """
    new_x_velocity = initial_x_velocity.*cos(angle_offset) - initial_z_velocity.*sin(angle_offset)
    new_z_velocity = initial_x_velocity.*sin(angle_offset) + initial_z_velocity.*cos(angle_offset)

    return (new_x_velocity,new_z_velocity)
end



end
