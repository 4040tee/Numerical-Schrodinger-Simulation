import numpy as np
import matplotlib.pyplot as plt

# Numerical Integration Functions
# -------------------------------

def spectral_radius(A):
    '''
    calculates the maximum absolute eigenvalue of given square matrix

    parameters:
    A: square matrix

    return:
    maximum absolute eigenvalue
    '''

    if A.shape[0] != A.shape[1]:

        print("matrix is not square")

    else:

        eig = np.linalg.eig(A)[0]
        max = np.abs(eig).max()

        return max


def make_tridiagonal(N, b, d, a):
    '''
    creates a tridiagonal (square) matrix given three values, 
    and a row/column size

    parameters:
    N: dimension of square matrix
    b: value on lower diagonal
    d: value on center diagonal
    a: value on upper diagonal

    return:
    square tridiagonal matrix of size A
    '''

    assert type(N) is int, "matrix dimension must be integer"
    assert N > 1, "matrix dimension must be greater than or equal to 1"

    # if N < 1:

    #     print("invalid size")
    #     return  N

    m = np.zeros((N,N), dtype=float)
    
    values = [[-1,b], [0,d], [1,a]]

    for value in values:

        index = value[0]

        if index >= 0:
            i = index

        else:
            i = (-index) * N

        m[:N-index].flat[i::N+1] = value[1]

    return m

# periodic boundary conditions taken from NM4P schro.py

def create_hamiltonian(N, h_bar, mass, h, V):
    '''
    Creates a tridiagonal hamiltonian operator matrix for the schrodinger equation. 
    Sets values according to the coefficient value in the equation, which incorporates natural units, mass, and grid size.
    Includes potential at different positions by adding potential inputs on the main diagonal.

    parameters:
    N: number of spatial grid points
    h_bar: natural units
    mass: mass
    h: grid size
    V: spatial potential array

    return:
    hamiltonian operator matrix
    '''

    assert type(N) is int, "matrix dimension must be integer"
    assert h_bar > 0 and mass > 0 and h > 0, "constant values must be non-zero and positive" 

    coeff = -h_bar**2/(2*mass*h**2)

    ham = make_tridiagonal(N,coeff, -2*coeff, coeff)

    # set potential inputs along main diagonal

    for i in range(1,N-1) :

        ham[i,i] = ham[i,i] + V[i]   
    
    # first and last rows for periodic boundary conditions

    ham[0,-1] = coeff;   ham[0,0] = -2*coeff;     ham[0,1] = coeff
    ham[-1,-2] = coeff;  ham[-1,-1] = -2*coeff;   ham[-1,0] = coeff

    return ham

# taken from NM4P schro.py

def create_initial_wavefunction(N, x,  k0, sigma0, x0, i_imag = 1j):
    '''
    creates a gaussian wavefunction to function as an initial condition for the schrodinger equation

    parameters:
    N: number of spatial grid points
    x: coordinates of spatial grid points
    k0: constant coefficient
    sigma0: constant coefficient
    x0: initial position
    i_imag: complex number representation

    return:
    initial gaussian wave function
    '''

    assert type(N) is int and N > 0, "vector dimension must be integer and non-zero"
    assert len(x) > 0, "spatial grid must be not be empty"

    Norm_psi = 1./(np.sqrt(sigma0*np.sqrt(np.pi)))   # Normalization

    psi = np.empty(N,dtype=complex)

    for i in range(N) :

        psi[i] = Norm_psi * np.exp(i_imag*k0*x[i]) * np.exp(-(x[i]-x0)**2/(2*sigma0**2))

    return psi

# overall structure taken from NM4P schro.py

def calculate_psi(N, max_iter, numerical_matrix, psi, tau):
    '''
    computes the time integration of the schrodinger equation over a given time interval
    utilising a matrix structured according to either the ftcs or crank-nicholson method.
    begins iteration from a given initial condition, expected input is a gaussian wave function.

    parameters:
    N: number of spatial grid points
    max_iter: number of time integrations to execute
    numerical_matrix: matrix which computes numerical solutions, structured according to either the ftcs or crank-nicholson method
    psi: initial condition wave function
    tau: time step size

    return:
    2-D array containing schrodinger equation results, 1-D array containing times corresponding to integration results
    '''

    assert type(N) is int and N > 0, "solution dimension must be a positive integer"
    assert max_iter > 0 and type(max_iter) is int, "number of iterations must be a positive integer"
    assert numerical_matrix.shape[0] == N and numerical_matrix.shape[1] == N, "numerical integration matrix must have dimensions N X N"
    assert tau != 0, "time step must be non-zero"

    plot_iter = max_iter/8                      
    
    iplot = 0
    t = np.empty(max_iter) 
    
    psis = np.full((N,max_iter), None, dtype=complex)
    psis[:,0] = psi[:]
    

    for iter in range(max_iter - 1) :

        # Compute new wave function using the Crank-Nicolson scheme

        psi = np.dot(numerical_matrix,psi) 
    
        iplot += 1
        t[iplot] = tau*(iter+1)
        psis[:, iplot] = psi[:]

    return psis, t

def sch_eqn(nspace, ntime, tau, method='ftcs', length=200, potential = [], wparam = [10, 0, 0.5]):
    '''
    computes the time integration of the schrodinger equation over a given time interval.
    begins by creating a hamiltonian matrix, then based on which method is selected 
    produces a matrix of numerical solution, either of the ftcs or crank-nicholson method.
    if ftcs method is selected, calculates spectral radius of matrix to determine if
    solution is expected to be stable, and terminates integration if solution is expected to be unstable.
    creates initial condition of time integration as a gaussian wave function according to given parameters.
    
    parameters:
    nspace: number of spatial grid points
    ntime: number of numerical time integrations to be computed
    tau: time step size, not relative time step
    method: method of solution, either 'ftcs', or 'crank', defaults to ftcs
    length: length of system, defaults to 200
    potential: array of spatial positions within the length of system which receive potential inputs
    wparam: parameters [sigma0, k0, x0] which function as coefficients in the gaussian equation of the initial condition
    -sigma0: standard deviation of the wavefunction
    -k0: average wavenumber
    -x0: location of the center of the wavepacket

    return:
    2-D array containing schrodinger equation results
    1-D array containing coordinates of spatial grid points
    1-D array containing times corresponding to integration results
    '''

    assert type(nspace) is int and nspace > 0, "number of spatial grid points must be a positive integer"
    assert ntime> 0 and type(ntime) is int, "number of iterations must be a positive integer"
    assert tau != 0, "time step must be non-zero"
    assert method == "ftcs" or method == "crank", "method must be either 'ftcs' or 'crank'"
    assert length > 0, "length must be non-zero and positive"

    # variables to be used in integration
    # -----------------------------------

    i_imag = 1j                 # imaginary unit
    N = nspace                  # number of spatial grid points
    L = float(length)           # length of system
    h = L/(N-1)                 # size of grid
    x = np.arange(N)*h - L/2.   # coordinates of spatial grid points
    h_bar = 1.                  # natural unit
    mass = 0.5                   # nass of particle
    tau = tau                   # size of time step

    velocity = 0.5              # velocity of particle

    max_iter = int(ntime)       # number of time integrations to be computed

    # create array V of positions of potential inputs

    V = np.zeros(N) # initialize to zeros

    # set positions where potential inputs are received to 1

    for v in potential:
        V[v] = 1
    
    # create hamiltonian operator matrix
    # ----------------------------------

    ham = create_hamiltonian(N, h_bar, mass, h, V)

    # based on selected method, compute numerical matrix
    # --------------------------------------------------

    method = method.lower()

    if method == "crank": # compute matrix for crank-nicholson method

        numerical_matrix = np.dot( np.linalg.inv(np.identity(N) + .5*i_imag*tau/h_bar*ham), 
                    (np.identity(N) - .5*i_imag*tau/h_bar*ham) )
    
    elif method == "ftcs": # compute matrix for ftcs method

        numerical_matrix = np.identity(N) - (i_imag*tau/h_bar)*ham

        # determine stability of solution, and end integration if solution is expected to be unstable

        if spectral_radius(numerical_matrix) > 1: 

            print("WARNING: solution is expected to be unstable. Returning zero results.")
            
            return np.zeros((N,max_iter+1)), x, np.zeros(max_iter+1)
        
    else: # end integration if an invalid method was entered as a parameter

        print("Invalid method, please enter either 'crank' or 'ftcs' as a method.")

        return np.zeros((N,max_iter+1)), x, np.zeros(max_iter+1)
    
    # return spectral_radius(numerical_matrix)

    # create the initial condition gaussian wave function 
    # ---------------------------------------------------

    # process the parameters entered by the user into their appropraite data types

    x0 = float(wparam[1])       # location of the center of the wavepacket
    k0 = float(wparam[2])       # average wavenumber
    sigma0 = float(wparam[0])   # standard deviation of the wavefunction
    
    # create the wave function

    psi = create_initial_wavefunction(N, x,  k0, sigma0, x0, i_imag)


    # compute the time integration over the interval given by the user
    # -----------------------------------------------------------------

    psis, t = calculate_psi(N, max_iter, numerical_matrix, psi, tau)
    
    return psis, x, t

# Plotting Functions
# -------------------

def plot_wave_func(psi, x , t, target_time = None):
    '''
    attempts to find an approximate match for a specified time 
    within the range of times from a numerically solved solution of the schrodinger equation, 
    and plots the wave function at that time.
    if specified time exceeds the range of times in the solution, plots the initial wave function.
    if no time is specified plots the initial wave function.
    
    parameters:
    psi:  2-D array containing schrodinger equation results
    x: 1-D array containing coordinates of spatial grid points
    t: 1-D array containing times corresponding to integration results
    target_time: time which should be attempted to plot

    output:
    plot of wave function
    '''

    assert len(psi) > 0 and len(x) > 0 and len(t) > 0, "solution must be non-empty"
 
    # if not time is specified set the time to 0

    if target_time is None:
        target_time = 0
    
    # initialize the time list index to 0

    time = 0

    # iterate through the range of times in the solution and
    # attempt to find an approximate match to the requested time

    try:

        for i in range(len(t)):

            if t[i] <= target_time and t[i+1] >= target_time:

                time = i
                break
        
        if i == 0:
            plt.title("Initial Wave Function")
        else:
            plt.title('Wave Function at time = ' + "{:.2f}".format(t[time]))

    # if the requested time exceeds the range of times, use the initial time

    except IndexError:

        print("wave function plot: given time exceeds the range of times")

        time = 0
        plt.title("Initial Wave Function")

    # plot the wave function

    plt.plot(x,np.real(psi[:,time])) # to add the imaginary part: ,'-',x,np.imag(psi[:,time]),'--')

    plt.xlabel('x')
    plt.ylabel(r'$\psi(x)$')
    plt.legend(('Real  ','Imag  '))
    plt.show()


def plot_prob_density(psi, x, t, target_time = None):
    '''
    attempts to find an approximate match for a specified time 
    within the range of times from a numerically solved solution of the schrodinger equation, 
    and plots the probability density at that time.
    probability density is calculated as the solution at that time multiplied by its conjugate
    if the specified time exceeds the range of times in the solution
    plots the initial probability density.
    if no time is entered plots all probability densities.

    parameters:
    psi:  2-D array containing schrodinger equation results
    x: 1-D array containing coordinates of spatial grid points
    t: 1-D array containing times corresponding to integration results
    target_time: time which should be attempted to plot

    output:
    plot of probability density
    '''

    assert len(psi) > 0 and len(x) > 0 and len(t) > 0, "solution must be non-empty"

    max_iter = np.shape(psi)[1] # define the range to iterate over as the length of the solutions

    psi_prob = np.real(psi * np.conjugate(psi)) # calculate the probability densities

    # initialize the time list index to 0

    time = 0 

    # if a time is specified, attempt to plot the density at that time
    
    if target_time is  not None:

        # iterate through the range of times in the solution and
        # attempt to find an approximate match to the requested time

        try:    
            for i in range(len(t)):
                
                print(t[i])

                if t[i] <= target_time and t[i+1] >= target_time:
                    
                    time = i
                    
                    break

            plt.title('Probability density at time = ' + "{:.2f}".format(t[time]))

         # if the requested time exceeds the range of times, use the initial time

        except IndexError:

            time = 0
            print("prob. density plot: given time exceeds the range of times")
            
            plt.title("Initial Probability Density")

        # plot the probability density

        plt.plot(x,psi_prob[:,time])

    # if a time is not specified plot all probability densities

    else:

        for col in range(max_iter):

            if psi[0,col] != None:

                plt.plot(x,psi_prob[:,col])

        plt.title('Probability density at various times')

    plt.xlabel('x')
    plt.ylabel('P(x,t)')
    plt.show()


def sch_plot(psi, x, t, plot="prob", target_time = None):
    '''
    given a numerical solution to the schrodinger equation,
    plots the probability density or wave function at the approximate match
    to the specified time.
    if time is not specified, or exceeds the range of times in the solution
    plots the initial density, or wave function.

    parameters:
    psi:  2-D array containing schrodinger equation results
    x: 1-D array containing coordinates of spatial grid points
    t: 1-D array containing times corresponding to integration results
    plot: type of plot to be produced, either 'prob' or 'psi'
    'prob' = probability density plot, 'psi' = wave function plot
    target_time: time which should be attempted to plot
    '''

    assert len(psi) > 0 and len(x) > 0 and len(t) > 0, "solution must be non-empty"
    assert plot == "prob" or plot == "psi", "plot method must be 'psi' or 'prob'"

    if plot == "psi":
        plot_wave_func(psi, x, t, target_time)
    elif plot == "prob":
        plot_prob_density(psi, x, t, target_time)
    else:
        print("Please enter either 'psi' for a wave function plot, or 'prob' for a probability density plot, as the 'plot' parameter.")

# Extra Functions
# ----------------

def get_spectral_radius_ftcs(nspace, tau, length=200, potential = []):
    '''
    creates a hamiltonian operator matrix, produces the numerical integration matrix for the ftcs solution,
    then calculates the spectral radius of the numerical integration matrix.
    if spectral radius exceeds 1, outputs a warning
    
    parameters:
    nspace: number of spatial grid points
    tau: time step size, not relative time step
    length: length of system, defaults to 200
    potential: array of spatial positions within the length of system which receive potential inputs
    wparam: parameters [sigma0, k0, x0] which function as coefficients in the gaussian equation of the initial condition
    -sigma0: standard deviation of the wavefunction
    -k0: average wavenumber
    -x0: location of the center of the wavepacket

    return:
    spectral radius of numerical integration matrix
    '''

    assert type(nspace) is int and nspace > 0, "number of spatial grid points must be a positive integer"
    assert tau != 0, "time step must be non-zero"
    assert length > 0, "length must be non-zero and positive"

    # variables to be used in integration
    # -----------------------------------

    i_imag = 1j                 # imaginary unit
    N = nspace                  # number of spatial grid points
    L = float(length)           # length of system
    h = L/(N-1)                 # size of grid
    x = np.arange(N)*h - L/2.   # coordinates of spatial grid points
    h_bar = 1.                  # natural unit
    mass = 0.5                   # nass of particle
    tau = tau                   # size of time step

    velocity = 0.5              # velocity of particle

    # create array V of positions of potential inputs

    V = np.zeros(N) # initialize to zeros

    # set positions where potential inputs are received to 1

    for v in potential:
        V[v] = 1
    
    # create hamiltonian operator matrix
    # ----------------------------------

    ham = create_hamiltonian(N, h_bar, mass, h, V)

    # based on selected method, compute numerical matrix
    # --------------------------------------------------

    numerical_matrix = np.identity(N) - (i_imag*tau/h_bar)*ham

    # determine stability of solution and return
    
    return spectral_radius(numerical_matrix)

def find_stable_tau(nspace, tau = 0.00001, length=200, potential = [], wparam = [10, 0, 0.5]):
    '''
    attempts to find a stable value of tau for the FTCS solution of the schrodinger equation,
    given inputs for an integration, and an initial value of tau.
    defaults tau value to 0.00001.
    
    parameters:
    nspace: number of spatial grid points
    tau: time step size, not relative time step. defaults to 0.00001
    length: length of system, defaults to 200
    potential: array of spatial positions within the length of system which receive potential inputs
    wparam: parameters [sigma0, k0, x0] which function as coefficients in the gaussian equation of the initial condition
    -sigma0: standard deviation of the wavefunction
    -k0: average wavenumber
    -x0: location of the center of the wavepacket

    return:
    stable value of tau
    '''

    assert type(nspace) is int and nspace > 0, "number of spatial grid points must be a positive integer"
    assert tau != 0, "time step must be non-zero"
    assert length > 0, "length must be non-zero and positive"

    print("Finding a stable tau value for tau...")

    spect = 2

    while spect > 1 and tau > 0:

        spect = get_spectral_radius_ftcs(nspace, tau, length, potential)

        if spect > 1: tau = tau / 2

    print("Stable tau value is: " + str(tau))

    return tau


