import numpy as np
import matplotlib.pyplot as plt

class VicsekModel:
    def __init__(self, N_particles, L, v, eta, r, dt):
        """
        Initialize the Vicsek model.
        
        Parameters:
        - N: Number of particles
        - L: Size of the square domain (L x L)
        - v: Speed of the particles
        - eta: Noise amplitude
        - r: Interaction radius
        - dt: Time step
        """
        self.N_particles = N_particles
        self.L = L
        self.v = v
        self.eta = eta
        self.r = r
        self.dt = dt
        
        # Initialize positions and orientations randomly
        self.positions = np.random.rand(N_particles, 2) * L
        self.orientations = np.random.rand(N_particles) * 2 * np.pi
        
        self.order_param = 0
    
    def update(self):
        """
        Update the positions and orientations of the particles.
        """
        # Compute new orientations
        new_orientations = np.zeros(self.N_particles)
        vel_directions = np.zeros( (self.N_particles,2) )
        
        summ = np.array([0.,0.])
        for i in range(self.N_particles):
            # Find neighbors
            dists = np.sqrt(np.sum((self.positions - self.positions[i])**2, axis=1))
            neighbors = np.where(dists < self.r)[0]
            
            # Compute average orientation of neighbors
            avg_orientation = np.arctan2(
                np.mean(np.sin(self.orientations[neighbors])),
                np.mean(np.cos(self.orientations[neighbors]))
            )
            
            # Add noise
            new_orientations[i] = avg_orientation + self.eta * (np.random.rand() - 0.5)
            
            vel_directions[i] = np.array( [np.cos(new_orientations[i]),\
                                     np.sin(new_orientations[i])] )
            summ += vel_directions[i]

        summ = np.linalg.norm(summ)
        self.order_param = summ/(self.N_particles)
        
        self.orientations = new_orientations

        # Update positions
        self.positions[:, 0] += self.v * np.cos(self.orientations) * self.dt
        self.positions[:, 1] += self.v * np.sin(self.orientations) * self.dt
        
        # Apply periodic boundary conditions
        self.positions = self.positions % self.L
    
    def simulate(self, steps):
        """
        Run the simulation for a given number of steps.
        
        Parameters:
        - steps: Number of simulation steps
        """
        for _ in range(steps):
            self.update()
         #   self.plot()
    
    def plot(self):
        """
        Plot the current state of the particles.
        """
        plt.figure(figsize=(8, 8))
        plt.quiver(
            self.positions[:, 0], self.positions[:, 1],
            np.cos(self.orientations), np.sin(self.orientations),
            angles='xy', scale_units='xy', scale=1
        )
        plt.xlim(0, self.L)
        plt.ylim(0, self.L)
        plt.show()
        
    def get_order_param(self):
        """
        This function gets the order parameter phi
        at a given time step.
        """
        return self.order_param
    
    def get_positions(self):
        return self.positions
    
    def get_velocities(self):
        """
        Returns: orientations
        -------
        TYPE: numpy array
            orientations of all particles

        """
        return self.orientations

# def main():
#     # Parameters
#     N_particles = 40        # Number of particles
#     L = 3.1        # Size of the domain
#     v = 0.03        # Speed of the particles
#     eta = 5     # Noise amplitude
#     r = 1.0         # Interaction radius
#     dt = 1.0        # Time step
#     steps = 100      # Number of simulation steps
    
#     # Initialize and run the Vicsek model
#     sim = VicsekModel(N_particles, L, v, eta, r, dt)
#     sim.simulate(steps)
    
#     phi = sim.get_order_param()
#     print(phi)
    #
#     r = sim.get_positions()
    
# if __name__ == "__main__":
#     main()
