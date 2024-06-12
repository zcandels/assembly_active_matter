import interacting_vicsek_particles as ap

domLen = 1
rho = 50
N_particles = int(rho*domLen**2)
dt = 1
T = 100
v = 0.03
interaction_radius= 1
eta = 0.1

Particles = ap.initialize_particles(v, N_particles)
    
ap.do_timestep(T, dt, domLen, Particles, interaction_radius, eta)

