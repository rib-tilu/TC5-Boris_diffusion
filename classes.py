# Vector and method classes for the implementation of the Boris algorithm

class Vector: ## for all 3D-vectors
    
    def __init__(self, x, y, z):
        
        self.x = x
        self.y = y
        self.z = z
        
    
    def __str__(self):
        
        return "(" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"
    
        
    def norm2(self):
        
        return (self.x**2 + self.y**2 + self.z**2)
    
    
        
    def norm(self):
        
        return self.norm2()**0.5
    
    
    
    def plus(self, other):
        
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
    
    
    
    def moins(self):
        
        return Vector(-self.x, -self.y, -self.z)
    
    
    
    def produit(self, other):
        
        return Vector(other*self.x, other*self.y, other*self.z)
    
    

## Functions used on objects from the Vector class ##

    

def produit_scalaire(u, v): # scalar product of 2 vectors
        
    return u.x*v.x + u.y*v.y + u.z*v.z



def produit_vectoriel(u, v): # cross product of 2 vectors
            
    return Vector(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x)
    


def gamma(v): # Lorentz factor calculus
    
    return 1/( 1 - produit_scalaire(v,v)/(299792458**2) )**0.5



def U(v): 
    
    return v.produit(gamma(v))



def gamma_inv(u):
    
    return 1/( 1 + produit_scalaire(u,u)/(299792458**2) )**0.5



def V(u):
    
    return u.produit(gamma_inv(u))



###############################################################################
###############################################################################

# Function taking as an input the position, velocity, fields vectors of the particle and performing a single 
# time step using the relativistic Boris algorithm, it returns the new position and velocity vectors


def step(posk, vitk, Bk, Ek, charge, masse, delta_t):
    
    pas = delta_t/2
    
    facteur = charge*pas/masse
    
    
    Um = U(vitk).plus(Ek.produit(facteur))
    
    T = (Bk.produit(facteur)).produit(gamma_inv(Um))
    
    inter = 2/(1 + produit_scalaire(T,T))
    
    
    S = T.produit(inter)
    
    Up = Um.plus(produit_vectoriel(Um.plus(produit_vectoriel(Um,T)), S))
    
    Uk1 = Up.plus(Ek.produit(facteur))
    
    
    Vk1 = V(Uk1)
    
    Xk1 = posk.plus(Vk1.produit(pas))
    
    return [Xk1, Vk1]