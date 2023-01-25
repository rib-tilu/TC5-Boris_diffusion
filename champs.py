# Fonctions donnant les expressions des champs E ou B utilisés

from classes import *
import numpy as np

###############################################################################
#                             CHAMP ÉLECTRIQUE                                #
###############################################################################

def E_const(E0):
    
    return Vector(-E0, 0, 0)


def E_noyau(Znoy, pos): # singe nucleus located at the origin
    
    return pos.produit((Znoy*1.602e-19)/(4*np.pi*8.85e-12*(pos.norm()**3)))


def E_plasma(Znoy, pos, d_i, l_d): # ion plasma Debye sphere field centered on the particle

    Ep = Vector(0,0,0)
    
    
    i_xl = int( (pos.x - l_d)//d_i )
    
    i_yl = int( (pos.y - l_d)//d_i )
    
    i_zl = int( (pos.z - l_d)//d_i )
    
    
    i_x = int( (pos.x + l_d)//d_i )
    
    i_y = int( (pos.y + l_d)//d_i )
    
    i_z = int( (pos.z + l_d)//d_i )
    
    
    
    for h in range(i_xl + 1, i_x + 1):
        
        for i in range(i_yl + 1, i_y + 1):
            
            for j in range(i_zl + 1, i_z + 1):
                
                
                ecart = pos.plus(Vector(h*d_i, i*d_i, j*d_i).moins())
                
                
                if ecart.norm() < l_d and h >=0:
                    
                    Ep = Ep.plus(E_noyau(Znoy, ecart))
    
    return Ep