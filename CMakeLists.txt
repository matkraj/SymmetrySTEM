find_package(ArrayFire)
add_executable(symmetry main.cpp symmetry_support.cpp)

# To use Unified backend, do the following.
# Unified backend lets you choose the backend at runtime
target_link_libraries(symmetry ArrayFire::afopencl)
