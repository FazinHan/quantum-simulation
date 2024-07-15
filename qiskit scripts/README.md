# quantum-simulation

First run `classical_solver.py`, the exact calculation. This saves an image to the `outputs` folder, and `qutip_data#.npz` to the `data` folder, where `#` is the latest number. 
To check: Does not run on GPU cluster.

The next script to be run is `vqe.py`, and this starts the VQE.

Finally, `main.py` plots the variational and exact calculations. Saves two images in `outputs` and one `dimer#.npz` to `data`.