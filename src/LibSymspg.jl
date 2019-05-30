module LibSymspg
using Libdl

# Load in `deps.jl`, complaining if it does not exist
const depsjl_path = joinpath(@__DIR__, "..", "deps", "deps.jl")
if !isfile(depsjl_path)
    error("LibSymspg not installed properly, run Pkg.build(\"LibSymspg\"), restart Julia and try again")
end
include(depsjl_path)

# Module initialization function
function __init__()
    # Always check your dependencies from `deps.jl`
    check_deps()
end

# Export our two super-useful functions
export call_libfoo

# Function to call into the `libfoo` shared library with the given arguments
function call_libfoo()
    # global libsymspg
    hdl = Libdl.dlopen_e(libsymspg)
    println(hdl)
    # hdl = Libdl.dlopen_e(libfoo)
    # @assert hdl != C_NULL "Could not open $libfoo"
    # foo = Libdl.dlsym_e(hdl, :foo)
    # @assert foo != C_NULL "Could not find foo() within $libfoo"
    # return ccall(foo, Cdouble, (Cdouble, Cdouble), a, b)
    return true
end

end #module LibSymspg
