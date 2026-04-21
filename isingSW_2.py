#!/usr/bin/env python
"""
Simulación del modelo de Ising 2D o 3D usando Swendsen-Wang (Optimizado con Numba)
"""
import platform
import argparse
import sys
import numpy as np
from numba import njit


@njit
def find(i, parent):
    root = i
    while parent[root] != root:
        root = parent[root]
    
    # Compresión de caminos
    curr = i
    while curr != root:
        nxt = parent[curr]
        parent[curr] = root
        curr = nxt
    return root

@njit
def union(i, j, parent):
    root_i = find(i, parent)
    root_j = find(j, parent)
    if root_i != root_j:
        parent[root_i] = root_j

def estado_inicial(L, dim, ordenado):
    if ordenado:
        if dim == 2:
            return np.ones((L, L), dtype=np.int8)
        else:
            return np.ones((L, L, L), dtype=np.int8)
    else:
        if dim == 2:
            return np.random.choice(np.array([-1, 1], dtype=np.int8), size=(L, L))
        else:
            return np.random.choice(np.array([-1, 1], dtype=np.int8), size=(L, L, L))

@njit
def energia_total(lattice, J, L):
    E = 0.0
    # Usamos .ndim en lugar de dim para que Numba compile bien
    if lattice.ndim == 2:
        for i in range(L):
            for j in range(L):
                spin = lattice[i, j]
                vecinos = lattice[(i+1)%L, j] + lattice[i, (j+1)%L] + lattice[(i-1)%L, j] + lattice[i, (j-1)%L]
                E += -J * spin * vecinos
    else:  # dim == 3
        for i in range(L):
            for j in range(L):
                for k in range(L):
                    spin = lattice[i, j, k]
                    vecinos = lattice[(i+1)%L, j, k] + lattice[(i-1)%L, j, k] + lattice[i, (j+1)%L, k] + lattice[i, (j-1)%L, k] + lattice[i, j, (k+1)%L] + lattice[i, j, (k-1)%L]
                    E += -J * spin * vecinos
    return E / 2.0
@njit
def magnetizacion(lattice):
    # Suma absoluta dividida por el tamaño
    return np.abs(np.sum(lattice)) / lattice.size


@njit
def swendsen_wang_simulacion(lattice, L, t_steps, T, J):
    p_enlace = 1.0 - np.exp(-2.0 * J / T)
    
    # Truco Numba: lattice.size nos da el total de espines sin tener que hacer L**dim
    N = lattice.size
    
    E_todas = np.zeros(t_steps + 1)
    M_todas = np.zeros(t_steps + 1)
    
    E_todas[0] = energia_total(lattice, J, L)
    M_todas[0] = magnetizacion(lattice)
    
    for t in range(t_steps):
        parent = np.arange(N)
            
        # CREAMOS LOS ENLACES
        if lattice.ndim == 2:
            for i in range(L):
                for j in range(L):
                    idx = i * L + j
                    if lattice[i, j] == lattice[i, (j+1)%L]:
                        if np.random.rand() < p_enlace:
                            union(idx, i * L + (j+1)%L, parent)
                    if lattice[i, j] == lattice[(i+1)%L, j]:
                        if np.random.rand() < p_enlace:
                            union(idx, ((i+1)%L) * L + j, parent)
        else: # dim == 3
            for i in range(L):
                for j in range(L):
                    for k in range(L):
                        idx = i*L*L + j*L + k
                        if lattice[i, j, k] == lattice[(i+1)%L, j, k]:
                            if np.random.rand() < p_enlace:
                                union(idx, ((i+1)%L)*L*L + j*L + k, parent)
                        if lattice[i, j, k] == lattice[i, (j+1)%L, k]:
                            if np.random.rand() < p_enlace:
                                union(idx, i*L*L + ((j+1)%L)*L + k, parent)
                        if lattice[i, j, k] == lattice[i, j, (k+1)%L]:
                            if np.random.rand() < p_enlace:
                                union(idx, i*L*L + j*L + ((k+1)%L), parent)

        # VOLTEAMOS CLÚSTERES
        cluster_spins = np.zeros(N, dtype=np.int8) 
        flat_lattice = lattice.flatten()
        
        for idx in range(N):
            root = find(idx, parent)
            if cluster_spins[root] == 0:
                cluster_spins[root] = 1 if np.random.rand() < 0.5 else -1
            
            flat_lattice[idx] = cluster_spins[root]
            
        lattice = flat_lattice.reshape(lattice.shape)

        E_todas[t+1] = energia_total(lattice, J, L)
        M_todas[t+1] = magnetizacion(lattice)
        
    return E_todas, M_todas

def main():
    parser = argparse.ArgumentParser(
        description='Simulación del modelo de Ising usando Swendsen-Wang',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('-L', '--L', type=int, required=True, help='Tamaño de la red (L x L)')
    parser.add_argument('-t', '--t_steps', type=int, required=True, help='Número de pasos de Swendsen-Wang')
    parser.add_argument('-Tmin', '--Tmin', type=float, required=True, help='Temperatura minima del sistema')
    parser.add_argument('-Tmax', '--Tmax', type=float, required=True, help='Temperatura máxima del barrido')
    parser.add_argument('-S', '-Step', type=int, required=True, help='Numero de elementos de temperatura')
    parser.add_argument('-J', '--J', type=float, required=False, help='Constante de acoplamiento', default=1.0)
    parser.add_argument('-d', '--d', type=int, required=True, help='Dimensión de la simulación')
    
    # CORRECCIÓN DE ARGPARSE
    parser.add_argument('--ord', action='store_true', help='Añade este flag para empezar ordenado')
    
    args = parser.parse_args()

    T_array = np.linspace(args.Tmin, args.Tmax, args.S)
    todos_los_datos = []
    
    # Imprimimos aviso
    print(f"Iniciando simulación Swendsen-Wang en {args.d}D para L={args.L} ({args.L**args.d} espines).")
    print(f"Estado inicial: {'Ordenado' if args.ord else 'Aleatorio'}\n")

    for T in T_array:
        # 1. Generamos red inicial
        lattice = estado_inicial(args.L, args.d, args.ord)
        
        # 2. Llamamos a la función JIT. (La primera vez tardará unos segundos extra en compilar, luego vuela)
        E_todas, M_todas = swendsen_wang_simulacion(
            lattice=lattice, L=args.L, t_steps=args.t_steps, T=T, J=args.J
        )

        pasos = np.arange(len(E_todas))
        T_columna = np.full(len(pasos), T)

        datos_esta_T = np.column_stack((T_columna, pasos, E_todas, M_todas))
        todos_los_datos.append(datos_esta_T)
        
        print(f"[OK] Simulada T = {T:.3f}")

    matriz_final = np.vstack(todos_los_datos)
    
    nombre_archivo = f"Datos_IsingSW_{args.d}D_L{args.L}_{args.ord}.txt"
    np.savetxt(nombre_archivo, matriz_final, 
               header="Temperatura Paso Energia_Total Magnetizacion_Espin", 
               fmt="%.4f %d %g %g", comments="# ")
    
    print(f"\n[¡COMPLETADO!] Archivo guardado: {nombre_archivo}")

if __name__ == '__main__':
    main()