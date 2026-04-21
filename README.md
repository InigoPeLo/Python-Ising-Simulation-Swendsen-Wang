# Swendsen-Wang Ising Model Simulator

A Python implementation of the 2D and 3D Ising model using the Swendsen-Wang cluster-update Monte Carlo algorithm. This approach effectively mitigates "critical slowing down" near phase transitions, offering faster thermal equilibration compared to standard local-update methods like Metropolis-Hastings.

* **Dimensionality:** Simulate both 2D square lattices and 3D cubic lattices.
* **Algorithm:** Implements the Swendsen-Wang cluster algorithm using a Disjoint Set Union (Union-Find) data structure to efficiently identify and flip clusters.
* **Visualization:** ```isingSW_ver2.py``` Includes options to display static 2D/3D snapshots of the lattice state at specific time steps, plus real-time animations tracking the lattice state, energy, and magnetization.
* **Data Export:** Automatically saves energy and magnetization time series to a text file for post-processing and analysis.

It also implements a configuration fix to ensure the animation works on Linux systems.

### Needed Arguments for ```isingSW_ver2.py```
| Argument | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `-L` | int | System size (Side Length) | Required |
| `-d` | int | System Dimension (2 or 3) | Required |
| `-t` | int | #steps of Montecarlo simulation | Required |
| `-T` | Float | System temperature | Required |
| `-J` | Float | Ferromagnetic interaction constant | 1.0 |
| `-p` | int | List of the time steps you want to plot (ex. `-p 0 50 100`) | 0 (None) |
| `-a` | Bool | Use `True` to see the evolution of the system animated | False |

### Needed Arguments for ```isingSW_2.py``` (numba optimized)
| Argument | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `-L` | int | System size (Side Length) | Required |
| `-d` | int | System Dimension (2 or 3) | Required |
| `-t` | int | #steps of Montecarlo simulation | Required |
| `-Tmin` | Float | Minimum temperature to start the thermal sweep | Required |
| `-Tmax` | Float | Maximum temperature for the sweep | Requiered |
| `-S` | int | Number of temperature steps between `Tmin` and `Tmax` | Required
| `-J` | Float | Ferromagnetic interaction constant | 1.0 |
| `-ord` | Bool | Initial state of the simulation True = cold start, False = hot start | False | 



For `isingSW_ver2.py` every time a simulation finishes, the script automatically generates a plain text file in the same folder. This file contains data recorded at each Monte Carlo step.

* **Filename format:** `isingSW-{Dimension}-{Lattice_Size}-{Temperature}.txt` 
  *(Example: `isingSW-2-50-2.269.txt`)*
* **Contents:** The file includes a header and two columns of data:
  1. **Energy ($E$)** of the system at step $t$.
  2. **Magnetization ($M$)** of the system at step $t$.

To run `ìsingSW_ver2.py`, follow these instructions:

1. **Clone the repository**
```bash
git clone https://github.com/InigoPeLo/Ising-Simulation-2D-3D-Swendsen-Wang-Animated.git
``` 
2. **Install requirements**
```bash
pip install -r requirements.txt
```
3. **Run the simulation**
You can run the simulation from the terminal adding the parameters of the simulation you want.

For example, for a system at $T=1.5$ with $L=8$ in 3D for a number of $50$ steps. Assuming you want to also run the simulation and that you dont need to plot any specific points.
```bash
python isingSW_ver2.py -L 8 -T 1.5 -d 3 -t 50 -a True -J 1
```
If on the other hand you want to plot the following time steps `0 10 25 50` you should run the following comand:
```bash
python isingSW_ver2.py -L 8 -T 1.5 -d 3 -t 50 -a True -J 1 -p 0 10 25 50
```
 
A file named `isingSW-3-8-1.5.txt` will be saved containing the energy $E$ and the magnetization $M$ for all time steps. 

`isingSW_ver2.py` is better for visualization but is much slower tham `ìsingSW_2.py`. I would recomend using the second script if you want to study the 3D ising model.

For `isingSW_2.py` every time a simulation finishes, the script automatically generates a plain text file in the same folder. This file contains data recorded at each Monte Carlo step.
* **Filename format:** `Datos_IsingSW_{Dimension}_{Lattice_Size}_{initial_state}.txt`
* ** Contente: ** The file includes a header and four columns of data:
  1. **Temperature**
  2. **Montecarlo Step**
  3. **Energy**
  4. **Absolute Magnetization |M|**
 
### Additional Requirements
To run `isingSW_2.py` , you will need to have the `numba` library installed in addition to `numpy`:
```bash
pip install numpy numba
```
#### Example of use
If you want to run a simulation with $L=32$ $t=10000$ $d=3$ and sweep between $T_{min}=2$ and $T_{max}=4$ with $S=10$ you need to run the following command:
```bash
python isingSW_2.py -L 32 -Tmin 2 -Tmax 4 -d 3 -t 10000 -S 10 -J 1
```
