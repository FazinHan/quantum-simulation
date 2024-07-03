def classical_plotter(num_qubits, B_arr, omega, energies, J, JII):
    for i in range(2**num_qubits):
        plt.plot(B_arr/omega, energies[:,i],'b.',ms=1.2)
    plt.xlabel('$B/\\Omega$')
    plt.ylabel('$\\varepsilon$')
    plt.title(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    

def qiskit_plotter(B_arr, singlets, triplets, omega, J, JII):
    # qiskit result plotter
    plt.plot(B_arr, singlets,'.')
    for i in range(3):
        plt.plot(B_arr, triplets[:,i],'.')
    plt.xlabel('$B/\\Omega$')
    plt.ylabel('$\\epsilon$')
    plt.title(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    

def qiskit_cost_plotter(B_arr, ls, costs, omega, J, JII):
    # costs plotter
    fig, axs = plt.subplots(np.unique(ls[:,0]).size,1) # number of layers is determined
    for idx, n in enumerate(ls[:,0]):
        locs = np.where(ls[:,0]==n)
        try:
            axs[idx].plot(ls[locs][:,1], costs[locs],'.')
            axs[idx].set_xlabel(f'layer {n[0]} steps')
            axs[idx].set_ylabel('cost')
        except TypeError:
            axs.plot(ls[locs][:,1], costs[locs],'.')
            axs.set_xlabel('steps')
            axs.set_ylabel('costs')
    fig.tight_layout()
    fig.suptitle(f'$\\Omega={omega}$, $J={J}$, $J_{{||}}={JII}$')
    