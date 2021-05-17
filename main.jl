include("./core/simulation.jl")

using CSV,DataFrames,Plots,HDF5
import Statistics

function april2021initials()
    """
    parameters i modify whenever for testing.
    """
    mass = 6.647e-27
    excitation_time = 1200e-9
    sep_zero = 0.0115  # y direction
    sep_gradient = (0.0297 - 0.0115)/0.105
    laser_zero =  0.0365 # x direction
    stray_zero = laser_zero
    x_velocity = 2000.0
    x_spread = 40.0
    z_velocity = 0.0
    z_spread = 0.0
    stray_constant = 0.0
    v_offset =-0.0
    angle_offset = 0.0#7.5*pi/180
    microwave_delay = 1000e-9

    return Initials(mass,x_velocity,x_spread,
    z_velocity,z_spread,
    excitation_time,
    angle_offset,
    laser_zero,sep_zero,sep_gradient,stray_zero,stray_constant,v_offset,microwave_delay)
end

function mycontrols()
    """
    The timing sequences of the world record run.
    """
    number_of_atoms = 20

    v_max = 1.0
    magnetic_field = 0.0e-4

    trace_file_path = "./data/inputs/april2021/test_longitudinal.csv"
    dipole_file_path = "./data/inputs/april2021/states_56_57.csv"

    return Controls(number_of_atoms,v_max,magnetic_field,trace_file_path,dipole_file_path)
end


function run_experiment_from_file()

    initials = april2021initials()
    controls = mycontrols()

    #volts = collect(range(-1.5,stop=1.5,length = 500))
    volts = collect(range(0.0,stop=5.0,length = 100))
    magnetic_field = 0.0
    prob = zeros(Float64,length(volts))
    stark = zeros(Float64,length(volts))
    dynamic = zeros(Float64,length(volts))
    seperation_phase = zeros(Float64,length(volts))
    seperation = zeros(Float64,length(volts))

    Base.Threads.@threads for i in 1:length(volts)
        initial_foo = deepcopy(initials)
        control_foo = deepcopy(controls)
        control_foo.v_max = volts[i]

        result = simulation(initial_foo,control_foo)
        prob[i] = result.probability_average
        stark[i] = result.stark_phase_average
        dynamic[i] = result.dynamic_phase_average
        seperation[i] = result.seperation_average
        seperation_phase[i] = result.seperation_phase_average

    end
    #p = plot(volts,0.5*(1.0.-cos.(stark+dynamic )),leg=false)
    p =plot(volts,dynamic,label="dynamic",legend=:bottomright)
    plot!(volts,stark,label="stark")
    plot!(volts,seperation_phase,label="seperation phase")
    plot!(volts,seperation_phase+dynamic+stark,label="total phase")
    plot!(volts,dynamic+stark,label="dyn + stark phase")

    #plot!(volts,dynamic+stark,label="total")
    g =plot(volts,seperation,label="sep",legend=:bottomright)


    #exp_data = CSV.read("src/data/march2021/data_1p000.csv",DataFrame)
    g=plot(volts,-1 .*prob .+1,label="total")
    plot!(volts,0.5.-0.5*cos.(seperation_phase+dynamic+stark),label="foo",legend=:bottomright,color="black")
    plot!(volts,0.5.-0.5*cos.(dynamic+stark),label="sep",legend=:bottomright)

    #plot!(volts,0.5*(1.0.-cos.(dynamic )),label="dynamic")

    display(p)
    display(g)
    return (p,g)
end
