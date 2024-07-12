import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from information import determine_next_filename
from physics import num_qubits as chain_length

'''$(i) A$(i_A)_length$(chain_length)_layer_plan$(layer_plan) $(timeevo) $(en)'''

# chain_length = 4

layer_plan = range(2,7) # if i is int, circuit will have i type 2 (complete) layers
layer_plan = [2] # if i is int, circuit will have i type 2 (complete) layers

# layer_plan = [4,5,[2,2,1,1,1],[2,1,1,1,2]]

# layer_plan = [1,2]
Jii = [.1,1]

if __name__=="__main__":
    B = float(sys.argv[1])
    def func(Jii):
        out = []
        for layers in layer_plan:
            if type(layers)==int:
                layers = [2]*layers
            pipe = os.popen(f'julia Julia_Algo1_noNoise.jl {chain_length} {B} {Jii} "{layers}"')
            out.append(pipe.read())
            # interest = out[1:-1]
            # print(len(out))
        return [float(i) for i in out]


    import time
    
    t0 = time.perf_counter()
    with ProcessPoolExecutor(len(Jii)) as exe:
        errors = [i for i in exe.map(func, Jii)]
    t1 = time.perf_counter()

    # print(errors)

    name = determine_next_filename(f'julia_result_errors_{B}_','txt','data')
    with open(name,'a+') as file:
        file.write(str(layer_plan)+'\n'+str(list(B))+'\n'+str(errors))
        print('written to',name)
    
    print(f'\ncomplete in {np.round(t1-t0,1)}s\n')
    os.system("git add .")
    os.system("git commit -m 'update'")
    os.system("git push")
    print('pushed successfully')
