# Swendsen-Wang Ising Model Simulator

A Python implementation of the 2D and 3D Ising model using the Swendsen-Wang cluster-update Monte Carlo algorithm. This approach effectively mitigates "critical slowing down" near phase transitions, offering faster thermal equilibration compared to standard local-update methods like Metropolis-Hastings.

* **Dimensionality:** Simulate both 2D square lattices and 3D cubic lattices.
* **Algorithm:** Implements the Swendsen-Wang cluster algorithm using a Disjoint Set Union (Union-Find) data structure to efficiently identify and flip clusters.
* **Visualization:** Includes options to display static 2D/3D snapshots of the lattice state at specific time steps, plus real-time animations tracking the lattice state, energy, and magnetization.
* **Data Export:** Automatically saves energy and magnetization time series to a text file for post-processing and analysis.

It also implements a configuration fix to ensure the animation works on Linux systems.

### Needed Arguments
| Argument | Type | Description | Default |
| :--- | :--- | :--- | :--- |
| `-L` | int | Tamaño del sistema (longitud de un lado de la red) | Requerido |
| `-d` | int | Dimensión del sistema (2 o 3) | Requerido |
| `-t` | int | Número total de pasos de Monte Carlo | Requerido |
| `-T` | Float | Temperatura del sistema | Requerido |
| `-J` | Float | Constante de acoplamiento ferromagnético | 1.0 |
| `-p` | int | Lista separada por espacios de los pasos para graficar estados estáticos (ej. `-p 0 50 100`) | 0 (Ninguno) |
| `-a` | Bool | Usa `True` para mostrar la animación en vivo | False |

Every time a simulation finishes, the script automatically generates a plain text file in the same folder. This file contains data recorded at each Monte Carlo step.

* **Filename format:** `isingSW-{Dimension}-{Lattice_Size}-{Temperature}.txt` 
  *(Example: `isingSW-2-50-2.269.txt`)*
* **Contents:** The file includes a header and two columns of data:
  1. **Energy ($E$)** of the system at step $t$.
  2. **Magnetization ($M$)** of the system at step $t$.

Tu run the script follow the next instructions

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
A file named `isingSW-3-8-1.5.txt` will be saved containing the energy $E$ and the magnetization $M$ for all time steps. 
