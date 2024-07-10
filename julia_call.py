import os
import matplotlib.pyplot as plt
import numpy as np

'''$(i) A$(i_A)_length$(chain_length)_layer_plan$(layer_plan) $(timeevo) $(en)'''

chain_length = 4

layer_plan = [4,5,[2,1,1,2]] # if i is int, circuit will have i type 2 (complete) layers

errors = []


for i in layer_plan:
    if type(i)==int:
        layers = [2]*i
    else:
        layers = i
    pipe = os.popen(f'julia Julia_Algo1_noNoise.jl {chain_length} 1 "{layers}"')
    out = pipe.readlines()[:]
    interest = out[1:-1]
    print(len(interest))
    errors.append(interest)

with open('julia_result.out','wb') as file:
    np.savez(file, layer_plan=layer_plan, errors=errors)

print('complete')