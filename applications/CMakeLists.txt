add_subdirectory(simall_beta)
add_subdirectory(simall_cavity)
add_subdirectory(simall_taylor_green)
add_subdirectory(simall_bfs)

# Week 20 — verification driver applications.  Each reads its reference
# values from tests/regression/cases/*.hpp and prints them, so operators
# can sanity-check the canonical acceptance numbers without launching the
# full GUI.  These executables ship in the runtime install component.
add_subdirectory(simall_naca)
add_subdirectory(simall_shocktube)
add_subdirectory(simall_rb)
add_subdirectory(simall_flameD)
add_subdirectory(simall_cylinder)
add_subdirectory(simall_pipe)
