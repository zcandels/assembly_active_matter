import unittest
#from vicsek_particles_flat import Particle

def PBC(self, p_x, p_y, L_x, L_y): 
    if (p_x < -L_x * 0.5):
        p_x = p_x + L_x
    if (p_x >= L_x * 0.5):
        p_x = p_x - L_x
    if (p_y < - L_y * 0.5):
        p_y = p_y + L_y
    if (p_y >= L_y * 0.5):
        p_y = p_y - L_y
    
    if abs(p_x) > L_x or abs(p_y) > L_y:
        raise Exception("Value out of bounds")
    return [p_x, p_y]

        

class TestPBC(unittest.TestCase):
    def test_outOfBounds(self):
        