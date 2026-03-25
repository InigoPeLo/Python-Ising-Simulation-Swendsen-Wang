#!/usr/bin/env python
"""
Simulación del modelo de Ising 2D o 3D usando Monte Carlo
"""
#Librerias
import matplotlib
import platform
import argparse
import sys

#Ventana externa en linux
if platform.system() == 'Linux':
    matplotlib.use('TkAgg')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches


#Definimos una clase para formar identificar clusters

class cluster:
    def __init__(self, size):
        #Al principio todos los spines forman un cluster
        self.parent = np.arange(size)

    def find(self, i):
        #Compresión de caminos, bucamos la raiz del cluster
        path = []
        while self.parent[i] != i:
            path.append(i)
            i = self.parent[i]
        for node in path:
            self.parent[node] = i
        return i
    #Unión de spines (creación de cluster)
    def union(self, i, j):
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            self.parent[root_i] = root_j




def estado_inicial(L, dim):
    """ Generación de estado inicial aleatorio  """
    if dim == 2:
        lattice = np.random.choice([-1, 1], size=(L, L))
    else:
        lattice= np.random.choice([-1, 1], size=(L, L, L))
    return lattice
    

# Función para calcular la energía total del sistema
def energia_total(lattice, J, L, dim):
    E=0
    if dim==2:
        for i in range(L):
            for j in range(L):
                spin= lattice[i, j]
                #Primeros vecinos con condiciones de frontera periódicas
                vecinos=lattice[(i+1)%L, j] + lattice[i, (j+1)%L] + lattice[(i-1)%L, j] + lattice[i, (j-1)%L]
                E+= -J * spin * vecinos
    else:  # dim == 3
        for i in range(L):
            for j in range(L):
                for k in range(L):
                    spin= lattice[i, j, k]
                    vecinos=lattice[(i+1)%L, j, k] + lattice[(i-1)%L, j, k] + lattice[i, (j+1)%L, k] + lattice[i, (j-1)%L, k] + lattice[i, j, (k+1)%L] + lattice[i, j, (k-1)%L]
                    E+= -J * spin * vecinos
    return E/2  # Cada par de interacciones se cuenta dos veces 

#Función calcular magnetización
def magnetizacion(lattice, dim):
    if dim==2:
        M = np.sum(lattice)/(len(lattice)**2)
    else:
        M=np.sum(lattice)/(len(lattice)**3)
    return M

#Definimos la función de Swendsen Wang

def swendsen_wang_ising(L, t_steps, T, dim, J=1, plot_indx=0):
#Con plot_indx escogemos que estados representar
    lattice= estado_inicial(L, dim)
    E= energia_total(lattice, J, L, dim)
    M= magnetizacion(lattice, dim)

    #Creamos listas de resultados (todos)
    E_todas= [E]
    M_todas= [M]
    lattice_todas= [lattice.copy()]

    #Creamos listas para animar
    E_anim=[E]
    M_anim=[M]
    lattice_anim=[lattice.copy()]
    
    #Probabilidad de formar un enlace
    p_enlace = 1.0 - np.exp(-2.0 * J / T)
    #Número total de espines
    N = L**dim
  
   
    for t in range(t_steps):
        cl = cluster(N)
            
        #Creamos los enlaces
        if dim == 2:
            for i in range(L):
                for j in range(L):
                    idx = i * L + j
                    #Solo comprobamos el vecino de abajo y el de la derecha para evitar repetir
                    #Vecino derecha
                    if lattice[i, j] == lattice[i, (j+1)%L]:
                        if np.random.rand() < p_enlace:
                            cl.union(idx, i * L + (j+1)%L)
                                
                        #Vecino abajo
                    if lattice[i, j] == lattice[(i+1)%L, j]:
                        if np.random.rand() < p_enlace:
                            cl.union(idx, ((i+1)%L) * L + j)
        #Dirección +1 en todos los ejes en 3d
        else:
            for i in range(L):
                for j in range(L):
                    for k in range(L):
                        idx = i*L*L + j*L + k
                        #+x
                        if lattice[i, j, k] == lattice[(i+1)%L, j, k]:
                            if np.random.rand() < p_enlace:
                                cl.union(idx, ((i+1)%L)*L*L + j*L + k)
                        #+y
                        if lattice[i, j, k] == lattice[i, (j+1)%L, k]:
                            if np.random.rand() < p_enlace:
                                cl.union(idx, i*L*L + ((j+1)%L)*L + k)
                        #+z
                        if lattice[i, j, k] == lattice[i, j, (k+1)%L]:
                            if np.random.rand() < p_enlace:
                                cl.union(idx, i*L*L + j*L + ((k+1)%L))

        cluster_spins = {}
        flat_lattice = lattice.flatten()
        
        for idx in range(N):
            root = cl.find(idx)
            #Si el clúster no tiene un espín asignado aún, se lo asignamos aleatoriamente +1 o -1
            if root not in cluster_spins:
                cluster_spins[root] = np.random.choice([-1, 1])
            flat_lattice[idx] = cluster_spins[root]
            
        lattice = flat_lattice.reshape(lattice.shape)

        #Recalculamos energía y magnetización
        E = energia_total(lattice, J, L, dim)
        E_t = E
        E_todas.append(E_t)
        M_todas.append(magnetizacion(lattice, dim))
        lattice_todas.append(lattice.copy())
        
        if t % max(1, (t_steps//100)) == 0:
            lattice_anim.append(lattice.copy())
            E_anim.append(E_t)
            M_anim.append(magnetizacion(lattice, dim))
            """Estado del sistema en diferentes momentos"""


    if plot_indx != 0:
            n_plots = len(plot_indx)
            filas = math.ceil(n_plots / 2)

            if dim == 2:
                fig, ax = plt.subplots(filas, 2, figsize=(10, 10))
                ax_flat = ax.flatten() if n_plots > 1 else [ax]
                fig.suptitle("Estado del sistema en momento t (Swendsen-Wang)")
                for i in range(len(plot_indx)):
                    im = ax_flat[i].imshow(lattice_todas[plot_indx[i]], cmap='RdBu', vmin=-1, vmax=1)
                    ax_flat[i].axis('off')
                    if plot_indx[i] == -1:
                        plot_indx[i] = t_steps
                    ax_flat[i].set_title(f"t = {plot_indx[i]}")

                    cbar = plt.colorbar(im, ax=ax_flat[i], ticks=[-1, 1], fraction=0.046, pad=0.04)
                    cbar.ax.set_yticklabels(['-1 ', '+1'])
                plt.show()

            else:
                fig, ax = plt.subplots(filas, 2, figsize=(10, 10), subplot_kw={'projection': '3d'})
                ax_flat = ax.flatten() if n_plots > 1 else [ax]
                fig.suptitle("Estado del sistema en momento t (3D Swendsen-Wang)")
                
                for i in range(len(plot_indx)):
                    matriz_3d = lattice_todas[plot_indx[i]]
                    colores = np.empty(matriz_3d.shape, dtype=object)
                    colores[matriz_3d == 1] = 'cyan'    
                    colores[matriz_3d == -1] = 'red'    
                    ax_flat[i].voxels(matriz_3d != 0, facecolors=colores, edgecolor='none', alpha=0.4)
                    ax_flat[i].axis('off')
                    
                    if plot_indx[i] == -1:
                        plot_indx[i] = t_steps
                        
                    ax_flat[i].set_title(f"t = {plot_indx[i]}")

                parch1 = mpatches.Patch(color='cyan', alpha=0.4, label='Espín +1')
                parch_men1 = mpatches.Patch(color='red', alpha=0.4, label='Espín -1')
                fig.legend(handles=[parch1, parch_men1], loc='upper right', title="Estados")
                plt.show()

    return lattice_todas, E_todas, M_todas, E_anim, M_anim, lattice_anim
    
def animacion(lattice_list, E_list, M_list, dim):

    if dim == 3:
        fig = plt.figure(figsize=(16, 5))
        ax_3d = fig.add_subplot(131, projection='3d')
        ax_e = fig.add_subplot(132)
        ax_m = fig.add_subplot(133)
        im = None
    else:
        fig, (ax_3d, ax_e, ax_m) = plt.subplots(1, 3, figsize=(15, 5))
        lattice_display = lattice_list[0]
        im = ax_3d.imshow(lattice_display, cmap='gray', vmin=-1, vmax=1)
        ax_3d.set_title('Estado de la red')
    
    # Para 3D, mostrar un cubo
    if dim == 3:
        L = lattice_list[0].shape[0]
        # Crear voxels para el estado inicial
        voxels_white = lattice_list[0] == 1
        voxels_black = lattice_list[0] == -1
        
        ax_3d.voxels(voxels_white, facecolors='white', edgecolor='k', linewidth=0.2)
        ax_3d.voxels(voxels_black, facecolors='black', edgecolor='k', linewidth=0.2)
        ax_3d.set_xlabel('X')
        ax_3d.set_ylabel('Y')
        ax_3d.set_zlabel('Z')
        ax_3d.set_title('Estado de la red (3D)')
        ax_3d.set_xlim(0, L)
        ax_3d.set_ylim(0, L)
        ax_3d.set_zlim(0, L)
    
    line_e, = ax_e.plot([], [], color='blue')
    ax_e.set_title('Energía')
    ax_e.set_xlim(0, len(E_list))
    ax_e.set_ylim(min(E_list), max(E_list))
    ax_e.set_xlabel('Pasos')
    ax_e.set_ylabel('E')
    
    line_m, = ax_m.plot([], [], color='red')
    ax_m.set_title('Magnetización')
    ax_m.set_xlim(0, len(M_list))
    ax_m.set_ylim(min(M_list) - 0.1, max(M_list) + 0.1)
    ax_m.set_xlabel('Pasos')
    ax_m.set_ylabel('M')

    def update(frame):
        artists = []
        
        if dim == 3:
            # Limpiar y redibujar el cubo
            ax_3d.clear()
            L = lattice_list[frame].shape[0]
            voxels_white = lattice_list[frame] == 1
            voxels_black = lattice_list[frame] == -1
            
            ax_3d.voxels(voxels_white, facecolors='white', edgecolor='k', linewidth=0.2)
            ax_3d.voxels(voxels_black, facecolors='black', edgecolor='k', linewidth=0.2)
            ax_3d.set_xlabel('X')
            ax_3d.set_ylabel('Y')
            ax_3d.set_zlabel('Z')
            ax_3d.set_title('Estado de la red (3D)')
            ax_3d.set_xlim(0, L)
            ax_3d.set_ylim(0, L)
            ax_3d.set_zlim(0, L)
        else:
            # Actualizar imagen 2D
            lattice_display = lattice_list[frame]
            im.set_data(lattice_display)
            artists.append(im)
        
        x_data = range(frame + 1)
        line_e.set_data(x_data, E_list[:frame + 1])
        line_m.set_data(x_data, M_list[:frame + 1])
        artists.extend([line_e, line_m])
        
        return artists

    ani = animation.FuncAnimation(fig, update, frames=len(lattice_list), 
                                  interval=50, blit=False)
    
    plt.tight_layout()
    plt.show()
def main():
    parser = argparse.ArgumentParser(
        description='Simulación del modelo de Ising usando Swendsen-Wang',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python ising_sw.py --L 40 --d 3 -t 100 -T 1.5 -J 1 --p 0 50 100 -a False
  python ising_sw.py -L 60 --d 2 -t 200 -T 0.5 -J 1 --p 0 -a True
        """
    )

    parser.add_argument('-L', '--L', type=int, required=True, help='Tamaño de la red (L x L)')
    parser.add_argument('-t', '--t_steps', type=int, required=True, help='Número de pasos de Swendsen-Wang')
    parser.add_argument('-T', '--T', type=float, required=True, help='Temperatura del sistema')
    parser.add_argument('-J', '--J', type=float, required=True, help='Constante de acoplamiento', default=1)
    parser.add_argument('-d', '--d', type=int, required=True, help='Dimensión de la simulación')
    parser.add_argument('-p', '--p', type=int, default=0, nargs='+', help='Indices a representar. Ej: -p 0 100, solo 0 es ninguno')
    parser.add_argument('-a', '--a', type=bool, default=False, required=False, help='Función para animar el algoritmo. True: Se ejecuta la animación, False: No se ejectua' )
    
    args = parser.parse_args()

    # Si pasaron solo un 0, lo guardamos como entero en lugar de lista para mantener tu lógica
    plot_indices = args.p if args.p != [0] else 0

    lattice_todas, E_todas, M_todas, E_anim, M_anim, lattice_anim = swendsen_wang_ising(
        L=args.L, dim=args.d, t_steps=args.t_steps, T=args.T, J=args.J, plot_indx=args.p
    )
    if args.a== True:
        animacion(lattice_anim, E_anim, M_anim, args.d)


 

    nombre_archivo = f"isingSW-{args.d}-{args.L}-{args.T}.txt"
    datos = np.column_stack((E_todas, M_todas))
    np.savetxt(nombre_archivo, datos, header="Energia Magnetizacion", fmt="%g", comments="# ")
    
    print(f"\n[ÉXITO] Resultados guardados en formato texto: {nombre_archivo}")


if __name__ == '__main__':
    main()

     