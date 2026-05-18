"""udf_velocity_profile.py — Parabolic inlet velocity UDF.

Demonstrates registering a vector-valued UDF that the solver can call at
every face centroid of the named inlet boundary.

Parabolic Poiseuille profile over a duct of half-height H, peak velocity
U_max along x.  y is the wall-normal coordinate.
"""
import simall

H      = 0.05      # m
U_max  = 2.0       # m/s

def inlet_velocity(t, x, y, z):
    # Time-independent parabolic profile.  Outside the duct -> zero.
    if abs(y) > H:
        return (0.0, 0.0, 0.0)
    u = U_max * (1.0 - (y / H) ** 2)
    return (u, 0.0, 0.0)

simall.UdfHost.instance().register_vector_txyz("inlet_velocity", inlet_velocity)
simall.log("Registered UDF 'inlet_velocity' (parabolic, U_max=%.2f m/s)" % U_max)
