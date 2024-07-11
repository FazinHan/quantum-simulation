include("Julia_Algo1_definitions.jl")

# use printed output to form plot: error = 1 + Overlap w/ time Evo

import SciPy as SP
import JLD2 as D2

import YaoPlots

#Parse the input
#setup for parallel -a Args.txt julia Julia_Algo1_noNoise.jl

chain_length = parse(Int64, ARGS[1])
# depth = parse(Int64, String(split(ARGS[1])[2]))
i_A = parse(Int64, ARGS[2])
JII = parse(Int64, ARGS[3]);

# chain_length = 2
# layers_string = rsplit(ARGS[3], "")
layer_plan = let expr = Meta.parse(ARGS[4])
    @assert expr.head == :vect
    Int.(expr.args)
end

J = 1.0;
omega = 5.0 * chain_length;
T = 2*π / omega;
time_iterations = 1000;
delta_t = T/(time_iterations);
times = LinRange(0,T,time_iterations+1);
times_trapezoid = LinRange(0+delta_t/2,T-delta_t/2,time_iterations);
A_vec = LinRange(0,1,10);
A = A_vec[i_A]

# for A in A_vec
println("For A = $(A)")

#Compute the zero state projector
psi_0 = Yao.zero_state(chain_length);
zero_state_projector =  Yao.matblock(Yao.state(psi_0)*transpose(Yao.state(psi_0)));

Num_EV = 2^chain_length <= 10 ? 2^chain_length : 10

EVs = zeros(Num_EV)

FausewehZhuCirc = getFausewehZhuCircuit(chain_length,layer_plan, U_T(chain_length,A, omega, delta_t, J, JII, times_trapezoid))

# println("circuit:")
YaoPlots.vizcircuit(FausewehZhuCirc, filename=joinpath(@__DIR__, "circ.png"))
   
number_of_parameters = Int((FausewehZhuCirc |> Yao.nparameters) / 2)
    
prev_solutions = []
    
Overlap_circ = getOverlapCircuit(chain_length,layer_plan)

en = 0

for i in 1:Num_EV

    optns = Dict("maxiter"=>6000, "disp"=>false)

    sol = nothing

    global en

    while true

        Theta = 2*pi*rand(number_of_parameters)

        cur_it = 0

        # function callbackfunc(xk)
        #     cur_it += 1
        #     if cur_it % 100 == 1
        #         println("Callback: $(VarL(xk, prev_solutions,chain_length, FausewehZhuCirc, Overlap_circ, zero_state_projector))")
        #     end
        # end
        #callbackfunc(xk) = return
        
        sol = SP.optimize.minimize(VarL, Theta, jac=dVarL, tol=1e-7, args=(prev_solutions,chain_length,FausewehZhuCirc,Overlap_circ, zero_state_projector), options=optns)

        (sol["success"] == true ) && break
    end
    
    #sol = SP.optimize.minimize(VarL, Theta,  args=(prev_solutions,chain_length,depth,A,FausewehZhuCirc,Overlap_circ), options=optns)

    en = en + Energy(sol["x"], chain_length, layer_plan, U_T(chain_length,A, omega, delta_t, J, JII, times_trapezoid), T)
    global timeevo = VarL(sol["x"], prev_solutions,chain_length, FausewehZhuCirc, Overlap_circ, zero_state_projector)
end

en = en/Num_EV;
    
println("A$(i_A)_length$(chain_length)_layers$(layer_plan) $(timeevo) $(en)")
flush(stdout)

append!(prev_solutions, [sol["x"]])
        
# end

println("Saving A$(i_A)_length$(chain_length)_depth$(layer_plan) ... ")
flush(stdout)

D2.save_object(".//results//Solution_A$(i_A)_length$(chain_length)_layers$(layer_plan).jld2", prev_solutions)



# end

# A = A_vec[1]

