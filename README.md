# quantum-simulation

First run `classical_main.py`, the exact calculation. This saves an image to the `outputs` folder, and `qutip_data#.npz` to the `data` folder, where `#` is the latest number. 
To check: Does not run on GPU cluster.

The next script to be run is `main.py`, and this either starts the VQE, or plots the variational and exact calculations (or both). Which one is run depends on which of the two ```if __name__="__main__":``` is chosen to be run. Saves two images in `outputs` and one `dimer#.npz` to `data`.