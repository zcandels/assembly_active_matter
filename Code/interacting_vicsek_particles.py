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
        
    def get_position(self, n):
        return self.position[n]
    
    def get_velocity(self, n):
        return self.velocity[n]
    
    def get_angle(self, n):
        return self.angle[n]
    
    def get_num_steps(self):
        position = self.position
        return len(position)
            
        
def initialize_particles(v, L_x, Nparticles):
    Particles = {}
    for i in range(Nparticles):
        r0 = np.random.uniform(-L_x, L_x, 2)
        angle0 = np.random.uniform(0, 2*np.pi)
        v0 = v*np.array([ np.cos(angle0), np.sin(angle0) ])
        Particles[i] = Particle(r0, v0, angle0)
    return Particles          


def do_timestep(T, dt, L_x, L_y, v, Particles, interaction_radius, eta):
    # Time-stepping
    Nt = int(round(T/float(dt)))
    Np = len(Particles)
    t = 0

    avg_normed_vel = 0

    for n in range(Nt-1):
        t += dt
        
        sum_vels = np.array([0., 0.])
        for i in range(Np):
            p = Particles[i]
            
            p_pos = p.get_position(n)
            
            avg_angle_p = 0
            ctr = 0
            for j in range( Np ):
                if i == j:
                    continue
                q = Particles[j]
                q_pos = q.get_position(n)
                
                pq_dist = p_pos - q_pos
                [x_dist, y_dist] = pq_dist

                dx = x_dist - L_x * round(x_dist/L_x)
                dy = y_dist - L_y * round(y_dist/L_y)

                r_pq2 = dx*dx + dy*dy
                if r_pq2 < interaction_radius*interaction_radius:
                    #print("interaction ", ctr)
                    ctr += 1
                    avg_angle_p += q.get_angle(n)
            
            avg_angle_p += p.get_angle(n)
            avg_angle_p /= ctr+1
            d_angle = np.random.uniform(-eta/2, eta/2)
            angle_np1 = avg_angle_p + d_angle
            p.set_angle(angle_np1)
            
            vel_np1 = v*np.array( [np.cos(angle_np1), np.sin(angle_np1)] )
            #print(np.linalg.norm(vel_np1))
            p.set_velocity(vel_np1)
            sum_vels += vel_np1

            p_vel = p.get_velocity(n)
            d_pos = p_vel*dt
            p.set_position(p_pos + d_pos)
            
        avg_normed_vel = np.linalg.norm(sum_vels)/(Np*v)
        return avg_normed_vel
        
    
def run_simulator():
    L_x = 10
    L_y = 10
    N_particles = 400
    rho = N_particles/(L_x*L_y)
    dt = 1
    T = 10000
    v = 0.03
    interaction_radius= 1
    eta = 5

    Particles = initialize_particles(v, L_x, N_particles)
    
    v_a = do_timestep(T, dt, L_x, L_y, v,
                      Particles, interaction_radius, eta)
    print(v_a)

def main():
    run_simulator()
    
if __name__ == "__main__":
    main()
