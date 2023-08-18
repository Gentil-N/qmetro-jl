

using QuantumOptics
using PyPlot
using LinearAlgebra

println("import done")
     
# Parameters
N_cutoff=10 #photon cutoff
N_atoms=6 #number of atoms
gamma=1.0     #decay rate
Δ=0     #detuning
g=10gamma    #coupling to the cavity
R=9gamma    #pump rate
κ=40gamma;    #rate of loss of photons from the cavity

# Bases
b_fock=FockBasis(N_cutoff) 
b_spin=SpinBasis(N_atoms//2)

# Fundamental operators
a = destroy(b_fock)
at = create(b_fock)
n = number(b_fock)

println(n)

sm = sigmam(b_spin)
sp = sigmap(b_spin)
sz = sigmaz(b_spin);

# Jaynes-Cummings-Hamiltonian
H0 = Δ*n
Hint = g*(at ⊗ sm + a ⊗ sp)
H = H0 ⊗ sparse(one(b_spin)) + Hint; 

# Initial state
Ψ0 = fockstate(b_fock, 0) ⊗ spindown(b_spin)

# Time interval
T_end=1
dt=0.01
T = [0:dt:T_end;]

# Collapse operators and decay rates
J=[sparse(one(b_fock)) ⊗ sm, sparse(one(b_fock)) ⊗ sp, a ⊗ sparse(one(b_spin))]
rates = [gamma, R, κ]

# Time evolution according to a master equation
tout, ρt = timeevolution.master(T, Ψ0, H, J; rates=rates);



# Photon number and inversion
exp_n_master= real(expect(n ⊗ sparse(one(b_spin)), ρt))
exp_sz_master = real(expect(sparse(one(b_fock)) ⊗ sz, ρt));

# Plot results
figure(figsize=(9, 3))
subplot(121)
ylim([0, N_cutoff])
plot(T, exp_n_master);
xlabel(L"\gamma t")
ylabel(L"\langle n \rangle")

subplot(122)
ylim([-N_atoms, N_atoms])
plot(T, exp_sz_master);
xlabel(L"\gamma t")
ylabel(L"\langle S_z \rangle")

tight_layout()

savefig("./res.svg")
     
# Photon number distribution
figure(figsize=(9,3))
ρ_end=ptrace(ρt[end], 2)
N=[0:1:N_cutoff;]
p_n=real(diag(ρ_end.data))
bar(N, p_n, label="photon number distribution")
legend()
xlabel("n")
ylabel(L"P(n)");

savefig("./pnd.svg")

# Spectrum
ρ_bar= (a ⊗ sparse(one(b_spin)))*ρt[end] 
τ_end=10
τ=collect(range(0.0, stop=τ_end, length=2^12))
τout, ρ_bar_τ=timeevolution.master(τ, ρ_bar, H, J; rates=rates) 
g_τ=expect(at ⊗ one(b_spin), ρ_bar_τ)
ω, spec = timecorrelations.correlation2spectrum(τ, g_τ; normalize_spec=true);
plot(ω, spec)
xlabel(L"\omega/\gamma")
ylabel("Intensity")
ylim(-0.1, 1.1)
xlim(-100, 100)
grid();

savefig("./spec.svg")
     
xvec=[-5:0.1:5;];
yvec=[-5:0.1:5;];