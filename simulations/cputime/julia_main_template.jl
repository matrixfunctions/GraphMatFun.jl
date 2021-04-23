using LinearAlgebra,BenchmarkTools,Printf;
USING

# Check CPU-turbo state.
# from https://askubuntu.com/questions/619875/disabling-intel-turbo-boost-in-ubuntu
# Disable:
# $ sudo -s
# # echo "1" > /sys/devices/system/cpu/intel_pstate/no_turbo
# Enable:
# $ sudo -s
# # echo "0" > /sys/devices/system/cpu/intel_pstate/no_turbo
turbo_disabled=NaN
if (isfile("/sys/devices/system/cpu/intel_pstate/no_turbo"))
    if readlines("/sys/devices/system/cpu/intel_pstate/no_turbo")[1] == "1"
        turbo_disabled=true;
    end
    if readlines("/sys/devices/system/cpu/intel_pstate/no_turbo")[1] == "0"
        turbo_disabled=false;
    end
end


if (turbo_disabled == false)
    @warn("CPU turbo is enabled. Timing likely unreliable.")
elseif (turbo_disabled == true)
    println("CPU turbo is disabled. Good!")
else
    println("CPU turbo state unknown.")
end


col=COLUMN;


n=2000;
A0=triu(tril(ones(n,n),3),-3)*1.0 +1.0*I;
A0[5,3] += 0.0001; # Break symmetry to avoid special case code
if col==1
    A0=2.5*A0/8
else
    A0=5.5*A0/8
end

println("Matrix norm: $(opnorm(A0,1))");
if (!isdefined(BLAS,:get_config))
    println("BLAS vendor: ",BLAS.vendor())
else
    println("BLAS config: ",BLAS.get_config())
end





### START REPEATED CODE
# Code to run: NAME *******
print(@sprintf("%20.20s: ","NAME"));
A=deepcopy(A0);
INCLUDE
bb=@benchmark FUNCTION($A);
mm=median(bb.times)*1e-9;
println("$(mm) mem: $(bb.memory)");
print("                       ")
@show bb.times

### END REPEATED CODE


println("DONE!");
