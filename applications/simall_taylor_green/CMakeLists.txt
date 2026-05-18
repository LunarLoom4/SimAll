add_executable(simall_taylor_green main.cpp)
target_link_libraries(simall_taylor_green PRIVATE
    SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Materials
    SimAll::Solver SimAll::Parallel SimAll::Io)
target_compile_features(simall_taylor_green PRIVATE cxx_std_20)
