"""
code which performs a velocity-verlet algorithim and at each step updates the dynamic and stark
phases of an array of atoms.
"""

module InternalLoop
using CSV,Interpolations,DataFrames

states = CSV.read("./data/inputs/april2021/states_56_57.csv",DataFrame)

field_lower = states.n1_field
energy_lower = states.n1_energy

field_higher = states.n2_field
energy_higher = states.n2_energy

intf_energy_lower = interpolate((field_lower,), energy_lower, Gridded(Linear()))
intf_energy_higher = interpolate((field_higher,), energy_higher, Gridded(Linear()))

mutable struct Atoms
    """
    Contains the positions of all the atoms in the bunch.
    These quantities get updated each time step.

    Parameters
    ==========
    x_positions: Vector{Float64}
        x positions of all the atoms
    x_velocities: Vector{Float64}
        x velocities of all the atoms
    x_accelerations: Vector{Float64}
        x accelerations of all the atoms
    z_positions: Vector{Float64}
        the z position through the plates for each atom
    z_velocities: Vector{Float64}
        the z velocities of each atom
    hmw_phase: Vector{Float64}
        the HMW phase accumulated in the z direction
    dynamic_phase_x: Vector{Float64}
        the dynamic phase accumulated in the x direction
    stark_phase: Vector{Float64}
        the total stark phase
    lower_energy_state: Bool
        determines which of the states the atom is in.
    mass: Float64
        mass of the atoms
    """

    x_positions::Vector{Float64}
    x_velocities::Vector{Float64}
    x_accelerations::Vector{Float64}

    z_positions::Vector{Float64}
    z_velocities::Vector{Float64}

    dynamic_phase::Vector{Float64}
    stark_phase::Vector{Float64}
    mass::Float64

    energy_getter :: Int
end

mutable struct Plates
    """
    The geometry and voltage of the force generating plates.

    applied_voltage: Float64
        the voltage on the plates
    magnetic_field: Float64
        applied magnetic field
    seperation_zero: Float64
        initial seperation of the plates
    seperation_gradient: Float64
        plate seperation gradient
    stray_zero: Float64
        position of canceled stray field
    stray_constant: Float64
        gradient of stray field
    """
    voltage::Float64
    magnetic_field::Float64
    sep_zero::Float64
    sep_gradient::Float64
    stray_zero::Float64
    stray_constant::Float64
end

function get_forces(x_positions::Vector{Float64},z_velocities::Vector{Float64},energy_getter::Int,plates::Plates)
    """
    Takes the x position of an atom and returns the force acting along this direction.
    """
    dx = 1e-8

    if energy_getter == false
        energy_function = intf_energy_lower
    elseif energy_getter == true
        energy_function = intf_energy_higher
    end
    ## get the electric field before and after
    elecs_before = get_electric_field.(x_positions .- dx,z_velocities,Ref(plates))
    elecs_after = get_electric_field.(x_positions .+ dx,z_velocities,Ref(plates))

    ## get the energy before and after
    energys_before  = energy_function.(elecs_before)
    energys_after = energy_function.(elecs_after)

    return -1.0*(energys_after .- energys_before)./(2*dx)
end

function get_electric_field(x_position::Float64,z_velocity::Float64,plates::Plates)
    voltage = plates.voltage
    plate_constant = plates.sep_gradient
    plate_sep_zero = plates.sep_zero
    stray_field_zero = plates.stray_zero
    stray_field_gradient = plates.stray_constant

    electric_field = (voltage/(plate_sep_zero + plate_constant*x_position)) + stray_field_gradient*(x_position-stray_field_zero) + plates.magnetic_field*z_velocity

    return electric_field
end

function update_canonical(atoms::Atoms,plates::Plates,dt::Float64)
    """
    updates the canonical variables of the atoms. Performs a velocity verlet algorithim.
    """
    new_position = atoms.x_positions .+ atoms.x_velocities.*dt .+ atoms.x_accelerations.*dt.*dt.*0.5
    new_force = get_forces(new_position,atoms.z_velocities,atoms.energy_getter,plates)
    new_acceleration =  new_force./atoms.mass
    new_velocity = atoms.x_velocities .+ (atoms.x_accelerations .+new_acceleration).*0.5.*dt

    atoms.x_positions = new_position
    atoms.x_velocities = new_velocity
    atoms.x_accelerations = new_acceleration

    atoms.z_positions = atoms.z_velocities*dt

end

function update_phases(atoms::Atoms,plates::Plates,dt::Float64)

    atoms.dynamic_phase = atoms.dynamic_phase + 1.0*(0.5*(atoms.mass)/ 1.0545718176461565e-34)* atoms.x_velocities.*atoms.x_velocities.* dt

    if atoms.energy_getter == false
        energy_function = intf_energy_lower
    elseif atoms.energy_getter == true
        energy_function = intf_energy_higher
    end

    elec_fields = get_electric_field.(atoms.x_positions,atoms.z_velocities,Ref(plates))
    atoms.stark_phase = atoms.stark_phase + energy_function.(elec_fields).*(dt/1.0545718176461565e-34)
end

function update_step(atoms::Atoms,plates::Plates,dt::Float64)
    """
    takes the current state of the atoms and the plates and updates them
    """
    update_canonical(atoms,plates,dt)
    update_phases(atoms,plates,dt)
end

end
