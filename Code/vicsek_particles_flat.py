import numpy as np
#from rdfpy import rdf
#import matplotlib.pyplot as plt

# I think the difference between this version and the previous 
# is the simple omission of the if statement 
# determining if i == j. 

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
    
    def PBC(self, n, L_x, L_y):
        p_x = self.position[n][0]
        p_y = self.position[n][1]
        
        if np.linalg.norm([p_x, p_y], np.inf) > L_x/2:
            print("out of bounds prior to application of BCs \n")
            print("p_x = ", p_x, " p_y = ", p_y)
        
        if (p_x < -L_x * 0.5):
            p_x = p_x + L_x
        if (p_x >= L_x * 0.5):
            p_x = p_x - L_x
        if (p_y < - L_y * 0.5):
            p_y = p_y + L_y
        if (p_y >= L_y * 0.5):
            p_y = p_y - L_y
        
        new_position = np.array([p_x, p_y])
        if np.linalg.norm([p_x, p_y], np.inf) > L_x/2:
            print("Still out of bounds \n")
            print("p_x = ", p_x, " p_y = ", p_y)
        self.position[n] = new_position

    def get_current_position(self, n):
        return self.position[n]

    def get_current_velocity(self, n):
        return self.velocity[n]

    def get_current_angle(self, n):
        return self.angle[n]

    def get_num_steps(self):
        position = self.position
        return len(position)


def initialize_particles(v, L_x, Nparticles):
    Particles = {}
    for i in range(Nparticles):
        r0 = np.random.uniform(-L_x/2, L_x/2, 2)
        angle0 = np.random.uniform(0, 2*np.pi)
        v0 = v*np.array([ np.cos(angle0), np.sin(angle0) ])
        Particles[i] = Particle(r0, v0, angle0)
    return Particles          


L_x = 3.1
L_y = 3.1
N_particles = 40
rho = N_particles/(L_x*L_y)
dt = 1.0
T = 100
v = 0.03
interaction_radius= 1.0
eta = 5.

g_r_array = [] # Empty list to store radial distribution function
# at each time step. 
radii_array = []



Particles = initialize_particles(v, L_x, N_particles)


# Time-stepping
Nt = int(round(T/float(dt)))
Np = len(Particles)
t = 0
avg_vels = np.zeros(Nt)
avg_normed_vel = 0

for n in range(Nt-1):
    t += dt
    coords = np.zeros( (Np, 2) )

    sum_vels = np.array([0., 0.])
    for i in range(Np):
        p = Particles[i]

        p_pos = p.get_current_position(n)

        coords[i] = p_pos

        sin_avg = 0 
        cos_avg = 0
        avg_angle_p = 0
        ctr = 0
        for j in range( Np ):
            q = Particles[j]
            q_pos = q.get_current_position(n)

            pq_dist = p_pos - q_pos
            # print(pq_dist)
            [x_dist, y_dist] = pq_dist

            dx = x_dist - L_x * round(x_dist/L_x)
            dy = y_dist - L_y * round(y_dist/L_y)

            r_pq2 = dx*dx + dy*dy
            if r_pq2 < interaction_radius*interaction_radius:
                #print("interaction ", ctr)
                ctr += 1
                q_ang = q.get_current_angle(n)
                sin_avg += np.sin(q_ang)
                cos_avg += np.cos(q_ang)

        sin_avg /= ctr
        cos_avg /= ctr
        avg_angle_p = np.arctan(sin_avg/cos_avg)
        d_angle = np.random.uniform(-eta/2, eta/2)
        angle_np1 = avg_angle_p + d_angle
        p.set_angle(angle_np1)

        vel_np1 = v*np.array( [np.cos(angle_np1), np.sin(angle_np1)] )
       # print(np.linalg.norm(vel_np1))
        p.set_velocity(vel_np1)
        sum_vels = sum_vels + vel_np1

        p_vel = p.get_current_velocity(n)
        d_pos = p_vel*dt
        p.set_position(p_pos + d_pos)

        p.PBC(n, L_x, L_y)

    #print(sum_vels)
    #[g_r, radii] = rdf(coords, dr = 0.01, parallel=False)#, progress=True)
    #g_r_array.append(g_r)
    #radii_array.append(radii)
    avg_normed_vel = np.linalg.norm(sum_vels)/(Np*v)
    avg_vels[n] = avg_normed_vel


