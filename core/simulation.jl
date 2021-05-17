include("./single_loop.jl")
include("./atom_dist.jl")
using DataFrames

mutable struct Initials
    mass :: Float64
    x_velocity :: Float64
    x_spread :: Float64
    z_velocity::Float64
    z_spread::Float64
    excitation_time::Float64
    angle_offset::Float64
    laser_zero::Float64
    sep_zero::Float64
    sep_gradient::Float64
    stray_zero::Float64
    stray_constant::Float64
    v_offset::Float64
    microwave_delay::Float64
end

mutable struct Controls
    """
    variables which can change be modified. including the paths for the files
    containing the experimental control traces and the dipoles.
    """
    number_of_atoms :: Int
    magnetic_field :: Float64
    v_max :: Float64
    trace_file_path :: String
    dipole_file_path :: String
end

struct InputTraces
    time_trace::Array{Float64,1}
    voltage_trace::Array{Float64,1}
    magnetic_trace::Array{Float64,1}
    dipole_trace1::Array{Float64,1}
    dipole_trace2::Array{Float64,1}
    dt::Float64
end

function generate_experiment(controls :: Controls)

    traces = CSV.read(controls.trace_file_path,DataFrame)

    times = traces.times
    dt = times[2] - times[1]
    voltages = traces.voltages * controls.v_max
    magnetic_trace = traces.magnetic * controls.magnetic_field
    dipole1 = traces.dipole1
    dipole2 = traces.dipole2
    output = InputTraces(times,voltages,magnetic_trace,dipole1,dipole2,dt)
    return output
end

function simulation(initials::Initials,controls::Controls)
"""
Given all the variables needed to run the experiment returns a result
"""


initial_canonicals = AtomDist.initial_distribution(initials.z_velocity,initials.z_spread,
                                                    initials.x_velocity,initials.x_spread,
                                                    initials.excitation_time,initials.microwave_delay,
                                                    initials.angle_offset,
                                                    controls.number_of_atoms,initials.laser_zero)


x_velocities = initial_canonicals[1]
z_velocities = initial_canonicals[2]
x_positions = initial_canonicals[3]
z_positions = initial_canonicals[4]

x_accelerations=zeros(Float64,controls.number_of_atoms)
hmw_phase = zeros(Float64,controls.number_of_atoms)
dyanmical_phase = zeros(Float64,controls.number_of_atoms)
stark_phase = zeros(Float64,controls.number_of_atoms)

my_atoms = SingleLoop.InternalLoop.Atoms(x_positions,x_velocities,x_accelerations,z_positions,z_velocities,dyanmical_phase,stark_phase,initials.mass,0)

res = generate_experiment(controls)

my_plates = SingleLoop.InternalLoop.Plates(controls.v_max,controls.magnetic_field,initials.sep_zero,initials.sep_gradient,initials.stray_zero,initials.stray_constant)

results = SingleLoop.run_experiment(my_atoms,my_plates,res.time_trace,res.voltage_trace,res.magnetic_trace,res.dipole_trace1,res.dipole_trace2,res.dt)
return results

end
