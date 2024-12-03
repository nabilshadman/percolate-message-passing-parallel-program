# 1. Percolate Message Passing Parallel Program

Percolate **message passing** parallel program using two-dimensional decomposition implemented in the C programming language.

This project utilizes the following tech stack:
- **C Programming Language:** The core of the simulation and its functionality.
- **Message Passing Interface (MPI):** For communication between processes in the parallel program.
- **GNU Make:** For compiling the code.
- **ImageMagick:** For visualizing the output maps.

These tools and technologies work together to perform parallel computations and visualizations efficiently.

## 2. Header and Source Files

In this section, we describe the **header** and **source** files of the percolate code.

**percolate.h**  
This is the header file for the simulation. System size (i.e. L) and the number of processes (i.e. NPROC) are defined here. Also, the file contains the function prototypes for visualizing the map and for generating random numbers.

**percolate.c**  
This source file contains the main function of the simulation. Communications with other processes (e.g. halo swaps, reduction operations) are carried out by calling routines from the Message Passing Interface (MPI).

**percio.c**  
This source file contains functions to visualize the map in the simulation. One version of the function is for dynamic arrays, and another is for static arrays.

**unirand.c**  
This source file contains functions for generating random numbers.

## 3. Compilation

In this section, we discuss how to **compile** the percolate code on the [Cirrus](https://www.epcc.ed.ac.uk/hpc-services/cirrus) supercomputer.

**3.1 The default configuration parameters for the percolate simulation are:**  
- System size: L = 768
- Number of processes: NPROC = 16
- Density of rock: rho = 0.4040

**3.2 Load relevant modules first:**  
Type in the command line:  
```bash
module load intel-compilers-19
module load mpt
module load ImageMagick
```

**3.3 To compile the code with default configuration parameters:**  
(1) Build the executable. Type in the command line:  
```bash
make
```

**3.4 To compile the code with different values of L and/or NPROC:**  
(1) Open `percolate.h` with any supported editor (e.g. Emacs).

(2) If changing L, change `#define L 768` --> `#define L <YOUR_L_OF_CHOICE>`.

(3) If changing NPROC, change `#define NPROC 16` --> `#define NPROC <YOUR_NPROC_OF_CHOICE>`.

(4) Save and exit.

(5) Build the executable. Type in the command line:  
```bash
make
```

**3.5 Running the simulation with a different value of rho**  
Any different value of rho (between 0.0 and 1.0) can be taken as an optional runtime argument. See the **Execution** section for more information on this.

## 4. Execution

In this section, we discuss our current **status** of the code, and how to **run** the code on both the frontend (i.e., login) and backend (i.e., compute) nodes of Cirrus.

**4.1 Note on the current version of the code**  
The current version of the code throws a segmentation fault on both login and backend nodes. Normally, you can use the command below to fix this:  
```bash
ulimit -s unlimited
```

However, the current version of the code throws a segmentation fault after the percolate update loop is complete and just after the map visualization IO begins (i.e., mapwrite). The segmentation fault is observed on both frontend and backend nodes. Additionally, the result of the percolate algorithm is not reliable on the frontend nodes for some configurations of rho, NPROC, and/or L.

Although the segmentation fault does occur on the backend nodes, the percolate algorithm completes properly and results in the correct answer (i.e., whether percolation occurs or not). So, we **recommend** running the current version of the code on the backend nodes for reliable correctness testing (i.e., without the visualization) and performance testing (e.g., measuring average time per step of the update loop).

A sample .out and .err file from a sample run of the percolate code on 2 backend nodes on Cirrus are included for reference in the repository. The sample run used the configuration below.

- L = 768
- NPROC = 64
- rho = 0.4040

We are currently working on debugging the subtle segmentation fault that occurs in the code. More information about this is included in the report.

**4.2 To run the code on the login node:**  
Type in the command line:  
```bash
mpirun -n <number_of_processes> ./percolate <seed> [rho]
```

Note that `<number_of_processes>` **must** match the `NPROC` value defined in `percolate.h`.

An integer value for `<seed>` is required (e.g., 7777).

Optionally, you may provide a float value of `[rho]` (e.g., 0.4039) after the seed. The value of rho **must** be between 0.0 and 1.0.

**Here are some examples:**  
*For NPROC = 16, seed = 7777, rho = 0.4040:*  
```bash
mpirun -n 16 ./percolate 7777
```

*For NPROC = 16, seed = 7777, rho = 0.4039:*  
```bash
mpirun -n 16 ./percolate 7777 0.4039
```

**4.3 To run the code on the backend node(s):**  
(1) Open `cirrusmpi.job` in any supported text editor (e.g., Emacs).

(2) If running on 1 node (i.e., 1-36 processes), change:  
`#SBATCH --ntasks=16` --> `#SBATCH --ntasks=<number_of_processes>`.

Note `<number_of_processes>` **must** match `NPROC` defined in `percolate.h`.

Note that each node on Cirrus has 36 CPU cores. So, for example, if you want to run the code with 72 processes, you need to run on 2 nodes.

If running on 2 nodes (i.e., 37-72 processes), also change:  
`#SBATCH --nodes=1` --> `#SBATCH --nodes=2`.

If running beyond 2 nodes (i.e., 73 or more processes), also change:  
`#SBATCH --qos=short` --> `#SBATCH --qos=standard`  
`#SBATCH --nodes=1` --> `#SBATCH --nodes=<number_of_nodes>`.

(3) If changing seed and/or rho, change:  
`srun --unbuffered --cpu-bind=core ./percolate 7777` -->  
`srun --unbuffered --cpu-bind=core ./percolate <seed> [rho]`.

(4) Save and exit.

(5) Type in the command line:  
```bash
sbatch cirrusmpi.job
```

This command will send the percolate executable to the slurm job scheduler, where the code will be run as soon as resources are available. The command will return a `<job_id>`. When the code is run, the output will be printed to a file titled `percolate-<job_id>.out`, and any error will be printed to a file titled `percolate-<job_id>.err`.

**Here are some examples (i.e., just the changes in cirrusmpi.job):**

*For NPROC = 16, seed = 7777, rho = 0.4040:*  
```bash
srun --unbuffered --cpu-bind=core ./percolate 7777
```

*For NPROC = 32, seed = 5555, rho = 0.4039:*  
```bash
#SBATCH --ntasks=32
srun --unbuffered --cpu-bind=core ./percolate 5555 0.4039
```

*For NPROC = 64, seed = 7777, rho = 0.4040:*  
```bash
#SBATCH --nodes=2
#SBATCH --ntasks=64
srun --unbuffered --cpu-bind=core ./percolate 7777
```

*For NPROC = 128, seed = 5555, rho = 0.4039:*  
```bash
#SBATCH --qos=standard
#SBATCH --nodes=4
#SBATCH --ntasks=128
srun --unbuffered --cpu-bind=core ./percolate 5555 0.4039
```

**4.4 Additional information**

For more information on running codes on the Cirrus system, please visit this [link](https://docs.cirrus.ac.uk/).
