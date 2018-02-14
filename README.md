# Optimization Algorithms
This repository contains some optimization algorithms in **Python**, namely **Particle Swarm Optimization Algorithm**, **Differential Evolution**, **Self Adaptive Differential Evolution** and **Differential Evolution with Composite Trial Vector Generation Strategies and Control Parameters**. All algorithms generates **logs** and **statitical results**, as well as **convegence and diversity graphs** for each run and overall runs, which is defined by the user.

-------------------------------------------------------------------------------------------------------------------------------------------
## Required Libs
matplotlib

## Usage
1. import the module (e.g import de)
2. create the object (e.g. d = de())
3. override fitness function with the function you want to maximize or minimize (e.g. d.fitness = YourFitnessFunction) OBS: you have to pass as parameter a **list** of the parameters you want to optimize with **size of dim**, and Rastrigin Function is implemented for testing purposes
4. then run the function optimizer (e.g. d.diferentialEvolution(parameters)) OBS: for PSO use d.particleSwarmOptimization(parameters)

## Parameters
1. pop_size: population size
2. dim: dimensions of the individual (parameters to optimize)
3. bounds: a tuple containing the bounds of the parameters (e.g. ((-5.12),(5.12)) bounds for 2 dimensions rastrigin benchmark function)
4. max_iterations: stop criteria
5. runs: number of algorithm executions for robustness
6. maximize: True to maximize and False to minimize (default: True)

## Algorithm Specific Parameters
### PSO
1. vmax: bounds for particle's velocities based on domain's percentage (default: 0.1 (10%))
2. c1: cognitive component (default: 2)
3. c2: social component (default: 2)
4. inertia: inertia weight for global and local search control (default: 0.7)

### DE
1. weight_factor (default: 0.8)
2. crossover_rate (default: 0.9)

### saDE
1. p1: the probability of applying strategy "rand/l/bin" (default: 0.5)
2. p2: the probability of applying strategy "current to best/2/bin" (default: 0.5)
3. learningPeriod: the probability of applying those two strategies update (default: 50)
4. crPeriod: CR values update (default: 5)
5. crmUpdatePeriod: CRm update (default: 25)

### coDE
1. param_pool: pool of parameters (suggested: param_pool =[[1.0,0.1], [1.0,0.9], [0.8,0.2]])

## References
Self-adaptive Differential Evolution Algorithm for Numerical Optimization (2005)
http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1554904

Differential Evolution with Composite Trial Vector Generation Strategies and Control Parameters (2011)
http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5688232

Particle Swarm Optimization
http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=488968

Differential Evolution – A Simple and Efficient Heuristic for global Optimization over Continuous Spaces
https://link.springer.com/article/10.1023%2FA%3A1008202821328

# Any question?
## gplichoski@gmail.com

