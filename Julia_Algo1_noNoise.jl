include("Julia_Algo1_definitions.jl")

import SciPy as SP
import JLD2 as D2

import YaoPlots

#Parse the input
#setup for parallel -a Args.txt julia Julia_Algo1_noNoise.jl

chain_length = parse(Int64, String(split(ARGS[1])[1]))
depth = parse(Int64, String(split(ARGS[1])[2]))
i_A = parse(Int64, String(split(ARGS[1])[3]))

Num_shots = 8192;
J = 1.0;
JII = 1.0;
omega = 5.0 * chain_length;
T = 2*π / omega;
time_iterations = 1000;
delta_t = T/(time_iterations);
times = LinRange(0,T,time_iterations+1);
times_trapezoid = LinRange(0+delta_t/2,T-delta_t/2,time_iterations);
A_vec = LinRange(0,2,20);
A = A_vec[i_A]

#Compute the zero state projector
psi_0 = Yao.zero_state(chain_length);
zero_state_projector =  Yao.matblock(Yao.state(psi_0)*transpose(Yao.state(psi_0)));

Num_EV = 2^chain_length

EVs = zeros(Num_EV)

FausewehZhuCirc = getFausewehZhuCircuit(chain_length,depth, U_T(chain_length,A, omega, delta_t, J, JII, times_trapezoid))

   
number_of_parameters = Int((FausewehZhuCirc |> Yao.nparameters) / 2)
    
prev_solutions = []
    
Overlap_circ = getOverlapCircuit(chain_length,depth)
    
for i in 1:Num_EV

    optns = Dict("maxiter"=>6000, "disp"=>false)

    sol = nothing

    while true

        Theta = 2*pi*rand(number_of_parameters)

        cur_it = 0

        function callbackfunc(xk)
            cur_it += 1
            if cur_it % 100 == 1
                println("Callback: $(VarL(xk, prev_solutions,chain_length, FausewehZhuCirc, Overlap_circ, zero_state_projector))")
            end
        end
        #callbackfunc(xk) = return
        
        sol = SP.optimize.minimize(VarL, Theta, jac=dVarL, callback=callbackfunc, tol=1e-7, args=(prev_solutions,chain_length,FausewehZhuCirc,Overlap_circ, zero_state_projector), options=optns)

        (sol["success"] == true ) && break
    end
    
    #sol = SP.optimize.minimize(VarL, Theta,  args=(prev_solutions,chain_length,depth,A,FausewehZhuCirc,Overlap_circ), options=optns)
    
    println("EV $(i) A$(i_A)_length$(chain_length)_depth$(depth) Overlap w/ time Evo = $(VarL(sol["x"], prev_solutions,chain_length, FausewehZhuCirc, Overlap_circ, zero_state_projector)) Energy = $(Energy(sol["x"], chain_length, depth, U_T(chain_length,A, omega, delta_t, J, JII, times_trapezoid), T))")
    flush(stdout)

    append!(prev_solutions, [sol["x"]])
    
end

println("Saving A$(i_A)_length$(chain_length)_depth$(depth) ... ")
flush(stdout)

D2.save_object(".//results//Solution_A$(i_A)_length$(chain_length)_depth$(depth).jld2", prev_solutions)



# end

# A = A_vec[1]

