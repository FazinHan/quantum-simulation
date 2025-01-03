import os
import sys
import numpy as np
#from concurrent.futures import ProcessPoolExecutor
from information import determine_next_filename


# layer_plan = [2,2] # if i is int, circuit will have i type 2 (dense) layers
# layer_plan = [2,[2,1],3,[1,1,1],[1,2,1],[1,1,1,2],[1,1,1,1],4] # if i is int, circuit will have i type 2 (dense) layers
layer_plan = [2,3,4,5]

Jii = [.1,.5,1]

if __name__=="__main__":
    chain_length = int(sys.argv[1])
    B = float(sys.argv[2])
    def func(Jii):
        out = []
        for layers in layer_plan:
            if type(layers)==int:
                layers = [2]*layers
            pipe = os.popen(f'/scratch/fizaank.phy21.iitbhu/julia-1.11.2/bin/julia Julia_Algo1_noNoise.jl {chain_length} {B} {Jii} "{layers}"')
            out.append(pipe.read().split())
            # interest = out[1:-1]
            # print(len(out))
        return [[float(i), float(j)] for i, j in out]


    import time
    
    t0 = time.perf_counter()
    data = np.array([func(i) for i in Jii])
    errors = data[:,0]
    grads = data[:,1]
    t1 = time.perf_counter()

    # print(errors)

    name = determine_next_filename(f'julia_result_errors_{chain_length}_{B}_','txt','data')
    with open(name,'a+') as file:
        file.write(str(layer_plan)+'\n'+str(list(Jii))+'\n'+str(errors))
        print('\nwritten to',name)
        name = determine_next_filename(f'julia_result_grads_{chain_length}_{B}_','txt','data')
    with open(name,'a+') as file:
        file.write(str(layer_plan)+'\n'+str(list(Jii))+'\n'+str(grads))
        print('\nwritten to',name)
    
    os.system(f'python jplot.py {chain_length}')
    print(f'\nplotted length {chain_length} and field {B}')
    
    print(f'\ncomplete in {np.round(t1-t0,1)}s\n')
    os.system("git add .")
    os.system("git commit -m 'update'")
    os.system("git push")
    print('pushed successfully')
