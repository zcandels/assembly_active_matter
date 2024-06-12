import numpy as np
import sys
import os
import time
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


class Particle():
    def __init__(self, r0, v0, angle0):
        self.position = []
        self.velocity = []
        self.angle = []
        
        self.position.append(r0)
        self.velocity.append(v0)
        self.angle.append(angle0)
        self.recent_side_collision = 0
        self.recent_top_bottom_collision = 0
        self.recent_corner_collision = 0

    
    def set_velocity(self,new_velocity):
        velocity = self.velocity
        velocity.append(new_velocity)
        return None

    def set_position(self, new_position):
        position = self.position
        position.append(new_position)
        return None
    
    def set_angle(self, new_angle):
        angle = self.angle
        angle.append(new_angle)
    
    def set_corner_collision_flag(self, n):
        self.recent_corner_collision_flag = n
        return None
        
    def set_side_collision_flag(self, n):
        self.recent_side_collision = n
        return None
    
    def set_topBottom_collision_flag(self, n):
        self.recent_top_bottom_collision = n
        return None
        
    def get_current_position(self, n):
        return self.position[n]
    
    def get_current_velocity(self, n):
        return self.velocity[n]
    
    def get_current_angle(self, n):
        return self.angle[n]
    
    def get_num_steps(self):
        position = self.position
        return len(position)
            
        
def initialize_particles(v, Nparticles):
    Particles = {}
    for i in range(Nparticles):
        r0 = np.random.uniform(-10, 10, 2)
        angle0 = np.random.uniform(0, 2*np.pi)
        v0 = v*np.array([ np.cos(angle0), np.sin(angle0) ])
        Particles[i] = Particle(r0, v0, angle0)
    return Particles          


def do_timestep(T, dt, v, Particles, interaction_radius, eta):
    # Time-stepping
    Nt = int(round(T/float(dt)))
    t = 0

    avg_normed_vel = 0

    for n in range(Nt-1):
        t += dt
        
        for i in Particles:
            p = Particles[i]
            
            p_pos = p.get_current_position(n)
            p_vel = p.get_current_velocity(n)
            
            d_pos = p_vel*dt
            
            p.set_position(p_pos + d_pos)
            
            v = np.linalg.norm(p_vel)
            
            avg_angle_p = 0
            ctr = 0
            for j in Particles:
                q = Particles[j]
                q_pos = q.get_current_position(n)
                
                dist = np.linalg.norm( p_pos - q_pos )
                if dist < interaction_radius:
                    ctr += 1
                    avg_angle_p += q.get_current_angle(n)
            
            avg_angle_p /= ctr+1
            d_angle = np.random.uniform(-eta/2, eta/2)
            angle_np1 = avg_angle_p + d_angle
            p.set_angle(angle_np1)
            
            vel_np1 = v*np.array( [np.cos(angle_np1), np.sin(angle_np1)] )
            p.set_velocity(vel_np1)
            avg_normed_vel += vel_np1
        avg_normed_vel = np.linalg.norm(avg_normed_vel)/len(Particles)
        
        
            

def output_results(Particles, domLen):
        [fig,ax] = plt.subplots()
        
    
def main():
    domLen = 1
    rho = 50
    N_particles = int(rho*domLen**2)
    dt = 1
    T = 100
    v = 0.03
    interaction_radius= 1
    eta = 0.1

    Particles = initialize_particles(v, N_particles)
    
    do_timestep(T, dt, domLen, Particles, interaction_radius, eta)
    
    output_results(Particles, domLen)
    
if __name__ == "__main__":
    main()
