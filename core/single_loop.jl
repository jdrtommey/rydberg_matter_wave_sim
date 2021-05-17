module SingleLoop
include("internal_loop.jl")

import Statistics

function run_arm(atoms::InternalLoop.Atoms,plates::InternalLoop.Plates,time_list::Array{Float64,1},voltage_list::Array{Float64,1},magnetic_list::Array{Float64,1},dipole_list::Array{Float64,1},dt::Float64)
    """
    given a voltage and dipole polarizability at each time step can compute the
    result for all atoms.
    """
    for (i,t) in enumerate(time_list)
        plates.voltage = voltage_list[i]
        plates.magnetic_field = magnetic_list[i]
        atoms.energy_getter = dipole_list[i]
        InternalLoop.update_step(atoms,plates,dt)
    end
end

struct RinterResult
    probability_average::Float64
    probability_std::Float64

    seperation_average::Float64
    seperation_std::Float64

    stark_phase_average::Float64
    stark_phase_std::Float64

    dynamic_phase_average::Float64
    dynamic_phase_std::Float64

    seperation_phase_average::Float64
end

function average_result(atoms::InternalLoop.Atoms,atoms2::InternalLoop.Atoms)

    seperation_phase = (atoms.mass/1.0545718176461565e-34)*(atoms2.x_velocities).*(atoms2.x_positions - atoms.x_positions)
    seperation_phase_average = Statistics.mean(seperation_phase)

    probability_average = Statistics.mean(0.5*(1.0.+cos.(seperation_phase.+(atoms.dynamic_phase-atoms.stark_phase) - (atoms2.dynamic_phase-atoms2.stark_phase))))
    probability_std = Statistics.std(0.5*(1.0.+cos.((atoms.dynamic_phase-atoms.stark_phase) - (atoms2.dynamic_phase-atoms2.stark_phase))))

    seperation_average = Statistics.mean(atoms.x_positions  - atoms2.x_positions )
    seperation_std = Statistics.std(atoms.x_positions  - atoms2.x_positions )

    stark_average = Statistics.mean(-atoms.stark_phase + atoms2.stark_phase)
    stark_std = Statistics.std(-atoms.stark_phase + atoms2.stark_phase)

    dynamic_average = Statistics.mean(atoms.dynamic_phase - atoms2.dynamic_phase)
    dynamic_std = Statistics.std(atoms.dynamic_phase - atoms2.dynamic_phase)

    return RinterResult(probability_average,probability_std,seperation_average,seperation_std,stark_average,stark_std,dynamic_average,dynamic_std,seperation_phase_average)
end

function run_experiment(atoms::InternalLoop.Atoms,plates::InternalLoop.Plates,time_list::Array{Float64,1},voltage_list::Array{Float64,1},magnetic_list::Array{Float64,1},dipole_list1::Array{Float64,1},dipole_list2::Array{Float64,1},dt::Float64)
    """
    Experiment measures the difference between two 'paths' in the interferometer, these are definded by
    the trace of the electric dipole moment as a function of time. Deep copies the initial state of the atoms
    and runs the experiment for the two different dipole time profiles. Returns a RinterResults struct
    which contains the averaged probabilities.
    """
    atoms2 = deepcopy(atoms)

    run_arm(atoms,plates,time_list,voltage_list,magnetic_list,dipole_list1,dt)
    run_arm(atoms2,plates,time_list,voltage_list,magnetic_list,dipole_list2,dt)

    return average_result(atoms,atoms2)
end

end
