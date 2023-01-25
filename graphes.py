# Fonctions de tracés graphiques 2D ou 3D

import numpy as np
import matplotlib.pyplot as plt



def plot_2d_graphs(listeX, listeY, titre, titreX, titreY):
    
    plt.plot(np.array(listeX), np.array(listeY), 'x', color='green')
    
    plt.xlabel(titreX)
    plt.ylabel(titreY)
    
    plt.title(titre)
    
    return None



def plot_2d_trajs(listeX, listeY, titre, titreX, titreY):
    
    plt.plot(np.array(listeX), np.array(listeY), color='red', label=titre)
    
    plt.xlabel(titreX)
    plt.ylabel(titreY)
    
    plt.title(titre)
    
    return None



def plot_3d(listeX, listeY, listeZ, titre, titreX, titreY, titreZ):
    
    ax = plt.axes(projection='3d')
    
    ax.plot3D(np.array(listeX), np.array(listeY), np.array(listeZ), 'gray')
    
    plt.title(titre)
    
    ax.set_xlabel(titreX, fontsize=13)
    ax.set_ylabel(titreY, fontsize=13)
    ax.set_zlabel(titreZ, fontsize=13)
    
    plt.show()
    
    return None

def plot_traj_noyau_2D(listeX, listeY, titre, titreX, titreY, rayon):
    
    plt.plot(np.array(listeX), np.array(listeY), color='red', label=titre)
    
    plt.plot(rayon*np.cos(np.linspace(0, np.pi, 100)), rayon*np.sin(np.linspace(0, np.pi, 100)), 'black')
    
    plt.xlim(-5e-8, 5e-8)
    plt.ylim(0, 6e-8)
    
    plt.xlabel(titreX)
    plt.ylabel(titreY)
    
    plt.title(titre)
    
    plt.show()
    
    return None
