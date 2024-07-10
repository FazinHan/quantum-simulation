import os
import matplotlib.pyplot as plt
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from information import determine_next_filename

'''$(i) A$(i_A)_length$(chain_length)_layer_plan$(layer_plan) $(timeevo) $(en)'''

chain_length = 4

layer_plan = [4,5,[2,1,1,2]] # if i is int, circuit will have i type 2 (complete) layers
# layer_plan = [1,2]

errors = []

if __name__=="__main__":
    def func(i):
        if type(i)==int:
            layers = [2]*i
        else:
            layers = i
        pipe = os.popen(f'julia Julia_Algo1_noNoise.jl {chain_length} 1 "{layers}"')
        out = pipe.readlines()
        interest = out[1:-1]
        print(len(interest))
        errors.append(interest)
        return 0

    with ProcessPoolExecutor(3) as exe:
        a = [0 for _ in exe.map(func, layer_plan)]

    name = determine_next_filename('julia_result_errors','txt','data')
    with open(name,'w') as file:
        file.write(str(layer_plan))
        file.write(str(errors.to_list()))
    
    print('complete')
    os.system("git commit -am 'update'")
    os.system("git push")
    print('pushed successfully')
