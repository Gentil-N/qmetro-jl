
using QuantumOptics
using PyPlot
using LinearAlgebra

println("import done")
     
# Parameters
gamma_b = 1.0 #1.054571817e−34 / 2

# Bases
b_spin = SpinBasis(1//2)

# Fundamental operators
sx = sigmax(b_spin);
sy = sigmay(b_spin);
sz = sigmaz(b_spin);

# Jaynes-Cummings-Hamiltonian
# H0 = Δ*n
# Hint = g*(at ⊗ a + a ⊗ a)
# H = H0 ⊗ sparse(one(b_spin)) + Hint;
H = -gamma_b*sz

println(H)

# Initial state
Ψ0 = 1/sqrt(2) * (spindown(b_spin) + spinup(b_spin))

println(Ψ0)

# Time interval
T_end=10
dt=0.01
T = [0:dt:T_end;]

# Time evolution according to a master equation
tout, ρt = timeevolution.schroedinger(T, Ψ0, H);

println(ρt)

dens_mat = tensor(Ψ0, dagger(Ψ0))

# Photon number and inversion
exp_sx_master = real(expect(dens_mat, ρt));
#exp_sy_master = real(expect(sy, ρt));
#exp_sz_master = real(expect(sz, ρt));

# Plot results
figure()
plot(T, exp_sx_master);
#plot(T, exp_sy_master);
#plot(T, exp_sz_master);
xlabel(L"t")
ylabel(L"prob")

tight_layout()

savefig("./res.svg")