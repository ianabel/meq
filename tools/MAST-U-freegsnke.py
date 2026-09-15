import matplotlib.pyplot as plt
import numpy as np


# build machine
from freegsnke import build_machine

from freegsnke import __path__ as path_to_freegsnke


tokamak = build_machine.tokamak(
    active_coils_path=f"{path_to_freegsnke[0]}/../machine_configs/MAST-U/MAST-U_like_active_coils.pickle",
    passive_coils_path=f"{path_to_freegsnke[0]}/../machine_configs/MAST-U/MAST-U_like_passive_coils.pickle",
    limiter_path=f"{path_to_freegsnke[0]}/../machine_configs/MAST-U/MAST-U_like_limiter.pickle",
    wall_path=f"{path_to_freegsnke[0]}/../machine_configs/MAST-U/MAST-U_like_wall.pickle",
)


from freegsnke import equilibrium_update
from freegs4e.boundary import freeBoundaryHagenow


eq = equilibrium_update.Equilibrium(
    tokamak=tokamak,      # provide tokamak object
    Rmin=0.1, Rmax=2.0,   # radial range
    Zmin=-2.2, Zmax=2.2,  # vertical range
    nx=513,                # number of grid points in the radial direction (needs to be of the form (2**n + 1) with n being an integer)
    ny=513,               # number of grid points in the vertical direction (needs to be of the form (2**n + 1) with n being an integer)
    boundary=freeBoundaryHagenow,
    # psi=plasma_psi
)

# plot the machine
import matplotlib.pyplot as plt

fig1, ax1 = plt.subplots(1, 1, figsize=(4, 8), dpi=80)
plt.tight_layout()

tokamak.plot(axis=ax1, show=False)
ax1.plot(tokamak.limiter.R, tokamak.limiter.Z, color='k', linewidth=1.2, linestyle="--")
ax1.plot(tokamak.wall.R, tokamak.wall.Z, color='k', linewidth=1.2, linestyle="-")

ax1.grid(alpha=0.5)
ax1.set_aspect('equal')
ax1.set_xlim(0.1, 2.15)
ax1.set_ylim(-2.25, 2.25)
ax1.set_xlabel(r'Major radius, $R$ [m]')
ax1.set_ylabel(r'Height, $Z$ [m]')


# initialise the profiles
from freegsnke.jtor_update import ConstrainPaxisIp
profiles = ConstrainPaxisIp(
    eq=eq,        # equilibrium object
    paxis=8e3,    # profile object
    Ip=6e5,       # plasma current
    fvac=0.5,     # fvac = rB_{tor}
    alpha_m=1.8,  # profile function parameter
    alpha_n=1.2   # profile function parameter
)


from freegsnke import GSstaticsolver
GSStaticSolver = GSstaticsolver.NKGSsolver(
    eq=eq,    # eq object
    # seed=42,  # random seed
    )


from freegsnke import GSstaticsolver

GSStaticSolver = GSstaticsolver.NKGSsolver(
    eq,
)

eq.tokamak.set_coil_current('Solenoid', 5000)
eq.tokamak['Solenoid'].control = False  # ensures the current in the Solenoid is fixed

import numpy as np
from freegsnke.inverse import Inverse_optimizer

Rx = 0.6      # X-point radius
Zx = 1.1      # X-point height
Rout = 1.4    # outboard midplane radius
Rin = 0.34    # inboard midplane radius

# set desired null_points locations (this can include X-point and O-point locations)
null_points = [[Rx, Rx], [Zx, -Zx]]

# set desired isoflux constraints with format
# isoflux_set = [isoflux_0, isoflux_1 ... ]
# with each isoflux_i = [R_coords, Z_coords]
isoflux_set = np.array([
    [[Rx, Rx, Rin, Rout, 1.0, 1.0, .8,.8], [Zx, -Zx, 0.,0., 2.0, -2.0, 1.62, -1.62]]
    ])

# set the coil current limits (upper and lower)
# coil ordering in this case: PX,  D1,  D2,  D3,  Dp,  D5,  D6,  D7,  P4,  P5,  P6
coil_current_limits = [
    [5e3, 9e3, 9e3, 7e3, 7e3, 5e3, 4e3, 5e3, 0.0, 0.0, None],
    [-5e3, -9e3, -9e3, -7e3, -7e3, -5e3, -4e3, -5e3, -10e3, -10e3, None]
]

# instantiate the freegsnke constrain object
constrain = Inverse_optimizer(
    null_points=null_points,
    isoflux_set=isoflux_set,
    coil_current_limits=coil_current_limits
)

# if you find coil limits are being violated or an undesireable solution is being produced,
# you can try increasing the penalty factor for violating the coil limits
# (here we just set it to its default value of 1e5)
constrain.mu_coils = 1e5

# solve!
GSStaticSolver.solve(eq=eq,
                     profiles=profiles,
                     constrain=constrain,
                     target_relative_tolerance=1e-6,
                     target_relative_psit_update=1e-3,
                     verbose=True, # print output
                     l2_reg=np.array([1e-12]*10+[1e-6]),
                     )


fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 8), dpi=80)

ax1.grid(zorder=0, alpha=0.75)
ax1.set_aspect('equal')
eq.tokamak.plot(axis=ax1,show=False)                                                          # plots the active coils and passive structures
ax1.fill(tokamak.wall.R, tokamak.wall.Z, color='k', linewidth=1.2, facecolor='w', zorder=0)   # plots the limiter
ax1.set_xlim(0.1, 2.15)
ax1.set_ylim(-2.25, 2.25)

ax2.grid(zorder=0, alpha=0.75)
ax2.set_aspect('equal')
eq.tokamak.plot(axis=ax2,show=False)                                                          # plots the active coils and passive structures
ax2.fill(tokamak.wall.R, tokamak.wall.Z, color='k', linewidth=1.2, facecolor='w', zorder=0)   # plots the limiter
eq.plot(axis=ax2,show=False)                                                                  # plots the equilibrium
ax2.set_xlim(0.1, 2.15)
ax2.set_ylim(-2.25, 2.25)


ax3.grid(zorder=0, alpha=0.75)
ax3.set_aspect('equal')
eq.tokamak.plot(axis=ax3,show=False)                                                          # plots the active coils and passive structures
ax3.fill(tokamak.wall.R, tokamak.wall.Z, color='k', linewidth=1.2, facecolor='w', zorder=0)   # plots the limiter
eq.plot(axis=ax3,show=False)                                                                  # plots the equilibrium
constrain.plot(axis=ax3, show=True)                                                          # plots the constraints
ax3.set_xlim(0.1, 2.15)
ax3.set_ylim(-2.25, 2.25)

plt.tight_layout()
