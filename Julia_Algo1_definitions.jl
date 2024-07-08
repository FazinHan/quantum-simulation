using LinearAlgebra
using Random

function pauli_operator(idx, num_sites, pauli_matrix)

    @assert(0 < num_sites)
    @assert(0 < idx <= num_sites) # Observe the difference!! Python: assert(0 <= idx < num_sites)

    gate_list = reduce(vcat, [
                                [Matrix{ComplexF64}(1.0I(2)) for _ in 1:idx - 1],        # [np.eye(2) for _ in range(idx)]
                                [Matrix{ComplexF64}(pauli_matrix)],                      # [pauli_matrix]
                                [Matrix{ComplexF64}(1.0I(2)) for _ in idx + 1:num_sites] # [np.eye(2) for _ in range(idx + 1, num_sites)]
                             ]
                      )

    return reduce(Base.kron, gate_list)
end;

function BField(A, t, omega)
    return [A*cos(omega*t), A*sin(omega*t)]
end;

function H_0(num_sites, J, JII)

    sigma_x = [pauli_operator(i,num_sites, [[0, 1.0] [1.0,0]]) for i in 1:num_sites];
    sigma_y = [pauli_operator(i,num_sites, [[0, -1.0im] [1.0im,0]]) for i in 1:num_sites];
    sigma_z = [pauli_operator(i,num_sites, [[1.0, 0.0] [0.0,-1.0]]) for i in 1:num_sites];

    return sum( vcat(   [J/4*sigma_x[i]*sigma_x[i+1] for i in 1:2:num_sites-1], 
                        [J/4*sigma_y[i]*sigma_y[i+1] for i in 1:2:num_sites-1], 
                        [J/4*sigma_z[i]*sigma_z[i+1] for i in 1:2:num_sites-1],
                        [JII*sigma_x[i]*sigma_x[i+2] for i in 1:2:num_sites-2],
                        [JII*sigma_y[i]*sigma_y[i+2] for i in 1:2:num_sites-2],
                        [JII*sigma_z[i]*sigma_z[i+2] for i in 1:2:num_sites-2],
                        [JII*sigma_x[(i+1)]*sigma_x[(i+1)+2] for i in 1:2:num_sites-3],
                        [JII*sigma_y[(i+1)]*sigma_y[(i+1)+2] for i in 1:2:num_sites-3],
                        [JII*sigma_z[(i+1)]*sigma_z[(i+1)+2] for i in 1:2:num_sites-3]
            ))
end;

function H_1(num_sites, time, A, omega)

    sigma_x = [pauli_operator(i,num_sites, [[0, 1.0] [1.0,0]]) for i in 1:num_sites];
    sigma_y = [pauli_operator(i,num_sites, [[0, -1.0im] [1.0im,0]]) for i in 1:num_sites];

    Bx, By = BField(A, time, omega)
    return sum( [Bx * sigma_x[i] + By * sigma_y[i] for i in 1:num_sites] )
end;

function U_T(num_sites, A, omega, delta_t, J, JII, times_trapezoid)
    H0 = H_0(num_sites, J, JII)
    U0 = LinearAlgebra.exp(-1.0im*delta_t*H0)

    U = I(2^num_sites)

    for t in times_trapezoid
        U *= U0*LinearAlgebra.exp(-1.0im*delta_t*H_1(num_sites,t,A, omega))
    end

    return U
end;

import Yao
import Yao.EasyBuild

function inverse_cnot(nbit, i, j)
    return Yao.cnot(nbit, j, i)
end

function entangle_map(num_sites, layer_type)
    return [[i i+1] for i in 1:2:num_sites-1]
end

function pair_ring_new(num_sites)
    return circshift(Yao.EasyBuild.pair_ring(num_sites), 1)
end

#Function to get the circuit
function getFausewehZhuCircuit(num_sites,depth, UT)

    var_circ = Yao.EasyBuild.variational_circuit(Float64, num_sites, depth)

    number_of_parameters = var_circ |> Yao.nparameters

    var_circ = Yao.dispatch!(var_circ,[k for k in 1:number_of_parameters]);

    var_circ_dag = Yao.dispatch!(Yao.chain(num_sites, Yao.put( (([i for i in 1:num_sites]) |> reverse) => Yao.EasyBuild.variational_circuit(Float64, num_sites, depth, entangle_map(num_sites,1), pair_ring_new(num_sites), entangler=inverse_cnot))),[-k for k in 1:number_of_parameters |> reverse]);

    circ = Yao.chain(num_sites, Yao.put(([i for i in 1:num_sites])=>var_circ), Yao.put(([i for i in 1:num_sites])=>Yao.matblock(UT, tag="UT")), Yao.put(([i for i in 1:num_sites])=>var_circ_dag) )

    return circ 

end

function updateFausewehZhuCircuit(circ, params)

    number_of_parameters = (circ |> Yao.nparameters)/2

    new_params = copy(params)

    append!(new_params, [-p for p in params] |> reverse)

    Yao.dispatch!(+, circ,new_params);
end

function setFausewehZhuCircuit(circ, params)

    number_of_parameters = (circ |> Yao.nparameters)/2

    new_params = copy(params)

    append!(new_params, [-p for p in params] |> reverse)

    Yao.dispatch!(circ,new_params);
end

function calcFausewehZhuGradient(grad)
    
    number_of_parameters = Int(size(grad)[1]/2)

    grad_true = grad[1:number_of_parameters]

    grad_true -= [p for p in grad[number_of_parameters+1:number_of_parameters*2] |> reverse]
    
    return grad_true
end

function getOverlapCircuit(num_sites,depth)

    var_circ = Yao.EasyBuild.variational_circuit(Float64, num_sites, depth, entangle_map(num_sites,1))

    number_of_parameters = var_circ |> Yao.nparameters

    var_circ = Yao.dispatch!(var_circ,[p for p in 1:number_of_parameters]);

    var_circ_dag = Yao.dispatch!(Yao.chain(num_sites, Yao.put( (([i for i in 1:num_sites]) |> reverse) =>Yao.EasyBuild.variational_circuit(Float64, num_sites, depth, entangle_map(num_sites,1), pair_ring_new(num_sites), entangler=inverse_cnot))),[-p for p in 1:number_of_parameters |> reverse]);

    circ = Yao.chain(num_sites, Yao.put(([i for i in 1:num_sites])=>var_circ), Yao.put(([i for i in 1:num_sites])=>var_circ_dag) )

    return circ
    
end

function setOverlapCircuit(circ, theta, theta_prime)

    param_true = vcat(theta, [-p for p in (theta_prime |> reverse) ])

    Yao.dispatch!(circ,param_true);
end

function OverlapWithOldSolutions(Theta, ListOfOldSolutions, num_sites, OverlapCirc, zero_state_projector)

    Cur_val = 0.0

    for Theta_prime in ListOfOldSolutions
        setOverlapCircuit(OverlapCirc,Theta, Theta_prime)
        Cur_val += Yao.expect(zero_state_projector, Yao.zero_state(num_sites)=>OverlapCirc)
    end

    return Cur_val
end



function VarL(Theta, ListOfOldSolutions, num_sites, FausewehZhuCirc, OverlapCirc, zero_state_projector )
    setFausewehZhuCircuit(FausewehZhuCirc, Theta)
    Cur_val = Yao.expect(zero_state_projector, Yao.zero_state(num_sites)=>FausewehZhuCirc)
            
    for Theta_prime in ListOfOldSolutions
        setOverlapCircuit(OverlapCirc,Theta, Theta_prime)
        Cur_val -= 50*Yao.expect(zero_state_projector, Yao.zero_state(num_sites)=>OverlapCirc)
    end

    return -real(Cur_val)
    
end;

function dVarL(Theta, ListOfOldSolutions, num_sites, FausewehZhuCirc, OverlapCirc, zero_state_projector )
    setFausewehZhuCircuit(FausewehZhuCirc, Theta)
    _, grad = (Yao.expect'(zero_state_projector, Yao.zero_state(num_sites)=>FausewehZhuCirc))
    grad_true = calcFausewehZhuGradient(grad)

    for Theta_prime in ListOfOldSolutions
        setOverlapCircuit(OverlapCirc,Theta, Theta_prime)
        _, grad2 = Yao.expect'(zero_state_projector, Yao.zero_state(num_sites)=>OverlapCirc)
        number_of_parameters = Int(size(grad2)[1]/2)
        grad2_true = grad2[1:number_of_parameters]
        grad_true -= 50*grad2_true
    end

    return -real(grad_true)
    
end;

function Energy(Theta, num_sites, depth, UT, T)

    var_circ = Yao.EasyBuild.variational_circuit(Float64, num_sites, depth)

    var_circ = Yao.dispatch!(var_circ,Theta)

    Result1 = Yao.zero_state(num_sites) |> var_circ 

    var_circ2 = Yao.chain(num_sites, Yao.put(([i for i in 1:num_sites])=>var_circ), Yao.put(([i for i in 1:num_sites])=>Yao.matblock(UT, tag="UT")))

    Result2 = Yao.zero_state(num_sites) |> var_circ2 

    state1 =  Result1 |> Yao.state

    state2 =  Result2 |> Yao.state

    return -angle(dot(state1, state2))/T

end