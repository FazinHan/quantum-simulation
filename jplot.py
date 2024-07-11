import os
import matplotlib.pyplot as plt
from information import determine_next_filename

name = determine_next_filename('julia_result_errors','txt','data',exists=True)
with open(name,'r') as file:
    layer_plan = file.readline()
    # layer_plan = eval(file.readline())
    # errors = eval(file.readline())
    errors = file.readline()

print(errors)
print(layer_plan)