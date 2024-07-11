import os
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from information import determine_next_filename

'''$(i) A$(i_A)_length$(chain_length)_layer_plan$(layer_plan) $(timeevo) $(en)'''

chain_length = 4
JII = range(0,1,.1)

layer_plan = [2,[2,1]] # if i is int, circuit will have i type 2 (complete) layers
# layer_plan = [1,2]
Jii = 0
B_idx = 1

errors = []

if __name__=="__main__":
    def func(layer_plan):
        if type(layers)==int:
            layers = [2]*layers
        pipe = os.popen(f'julia Julia_Algo1_noNoise.jl {chain_length} {B_idx} {Jii} "{layers}"')
        out = pipe.readlines()
        interest = out[1:-1]
        print(len(interest))
        errors.append(interest)
        return 0

    import time
    
    t0 = time.perf_counter()
    with ProcessPoolExecutor(len(layer_plan)) as exe:
        a = [0 for _ in exe.map(func, layer_plan)]
    t1 = time.perf_counter()

    name = determine_next_filename('julia_result_errors','txt','data')
    with open(name,'w') as file:
        file.write(str(layer_plan)+'\n'+str(errors))
        print('written to',name)
    
    print(f'complete in {np.round(t1-t0,1)}s')
    os.system("git commit -am 'update'")
    os.system("git push")
    print('pushed successfully')
