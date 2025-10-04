#= ====     ***** tutorial to show how it works *****     ==== =#
using PyCall
using PyPlot
using LinearAlgebra
using Polynomials
using FFTW

t=0:0.1:3;
t2=-1:0.01:4;

fnts(t)=sinpi.(2t) .+ 0.01;
ts = fnts(t);
ts2= fnts(t2);

M = 7; N = 5;
TotN = N + M;
# N > rank(M) will be required to obtain stable solution 
#     to avoid singularity by degeneration

fs=22;
figure(tight_layout=true);
xlim(-0.05,3.05);
plot([-0.1,3.1],[0,0],label="", color="black");
# plot(t2,ts2, label="sin2πt+0.01", lw=2);
plot(t2,ts2, label="0.01+sin2πt", lw=2);
scatter(t[1:TotN], ts[1:TotN], label="$TotN samples for analysis", color="red", lw=2);
xlabel("time (s)",fontsize=fs); 
ylabel("amplitude (arb. unit)",fontsize=fs); 
xticks(fontsize=fs); yticks(fontsize=fs);
legend(fontsize=16, loc="upper right"); 
annotate("(a)",(2.7,0.9),fontsize=16,zorder=6);
show();

myeps = 1e-15;

# ----------------------------------- internal functions;
"""find extreme value in matrix A""";
extreme(Aᵢⱼ) =
    Aᵢⱼ |> maximum |> abs > Aᵢⱼ |> minimum |> abs ? maximum(Aᵢⱼ) : minimum(Aᵢⱼ);
"""find exterme value in vector x""";
absmax(x) = maximum(abs.(extrema(x)));
"""corresponding index of above vector x""";
absmaxidx(x) = filter(i -> abs(x[i]) == absmax(x), eachindex(x));

# ----------------------------- set autocorrelation eq. AX=B
A = zeros(M, M);
for m in 1:M, m′ in 1:M, n in 0:N - 1
    A[m,m′] += ts[TotN - m - n] * ts[TotN - m′ - n]
end
A

B = zeros(M);
for m′ in 1:M, n in 0:N - 1
    B[m′] += ts[TotN - 0 - n] * ts[TotN - m′ - n]
end
B'

# ------------------------------------ normalize matrix
scaler = extreme(A); 
A /= scaler;
B /= scaler;
A
B'

# ----------------- care for void records: safety reason
for i in 1:M
    if norm(A[i,:]) < myeps
        A[i,:] = zeros(M)
        B[i] = 0
        A[i,i] = 1
    end
end
A
B'

# ------------------------------------ pivotting
for i in 1:M - 1
    amidx = absmaxidx(A[i:M,i])[1] + i - 1
    A[amidx,:], A[i,:] = A[i,:], A[amidx,:]
    B[amidx], B[i] = B[i], B[amidx]
end
A
B'

# ----------------------------------- sweepout forward
for i in 1:M - 1
    if abs(A[i,i]) < myeps
        A[i,:] = zeros(M)
        B[i] = 0
        A[i,i] = 1
    end
    for j = i + 1:M
        mx = A[j,i] / A[i,i]
        A[j,:] .-= mx * A[i,:]
        B[j] -= mx * B[i]
    end
end
A
B'

# --------------------------------- sweepout backward
for i in M:-1:2
    if abs(A[i,i]) < myeps
        A[i,:] = zeros(M)
        B[i] = 0
        A[i,i] = 1
    end
    for j = 1:i - 1
        mx = A[j,i] / A[i,i]
        A[j,:] .-= mx * A[i,:]
        B[j] -= mx * B[i]
    end
    B[i] /= A[i,i]
    A[i,:] /= A[i,i]
end
B[1] /= A[1,1];
A[1,:] /= A[1,1];
A
B'

# ------------------------- remove null higher orders aₘ (m>M′)
M′ = M;
for i in 1:M
    if abs(B[i]) > myeps
        M′ = i
    end
end
M′

# ------------------------------ set prediction coeffs. & get modes
predcoeffs = ones(M′ + 1);
for i in 1:M′
    predcoeffs[i] = -B[M′ + 1 - i]
end
predcoeffs'
Polynomial(predcoeffs)
modes = roots(Polynomial(predcoeffs))

# ------------ remove exterme modes which corrrespond to noise floor
MaxDiffBetweenEdges = 100;
M′′ = 0;
for i in 1:M′
    growwidth = abs(modes[i])^TotN
    if maximum([growwidth, 1 / growwidth]) < MaxDiffBetweenEdges
        M′′ += 1
        modes[M′′] = modes[i]
    end
end
modes = modes[1:M′′]
M′′

# ---------- calc complex amps. at left bound solve AX=B
A = zeros(ComplexF64, M′′, M′′);
B = zeros(ComplexF64, M′′);

for i in 1:M′′, j in 1:M′′, n in 0:TotN - 1
    A[i,j]  += (modes[i] * modes[j])^n
end
A

for i in 1:M′′, n in 0:TotN - 1
    B[i] += ts[n + 1] * modes[i]^n
end
B

iCAmp = A \ B

# -------------------  prepare return values
iFrq = imag(log.(modes)) / 2π;
iAVR = real(log.(modes));
results = [];
for i in 1:M′′
    push!(results, (iFrq = iFrq[i], iAVR = iAVR[i], iCAmp = iCAmp[i]))
end
results
results=reverse(results);

# ------------------- plot spectrum of each mode
""" equation for Lorentz profile spectrum """;
LorentzProf(f,iFrq,iAVR,iCAmp) =
    abs(iCAmp) * √(iAVR^2 / ((2π * (abs(iFrq) - f))^2 + iAVR^2));

    figure(tight_layout=true);
    #--------------------- FFT for demo
    yFFT=abs.(fft(ts[1:TotN]))[1:7]/TotN;
    xFFT=0:1:6;
    bar(xFFT/1.2,yFFT,width=1/1.2,color="#FFFFFF",edgecolor="#000000",label="Fourier", lw=2);
    #---------------------------
    x = 0:1e-4:0.5;
    for i in 1:M′′
    #   y = LorentzProf.(x, iFrq[i],         iAVR[i],         iCAmp[i]        );
    #   y = LorentzProf.(x, results[i][1],   results[i][2],   results[i][3]   );
        y = LorentzProf.(x, results[i].iFrq, results[i].iAVR, results[i].iCAmp);
        plot(10x, y, label="mode no. $i", lw=2);
    end
    annotate("(b)",(4.5,0.1),fontsize=16);
    xlabel("frequency (Hz)", fontsize=fs); 
    ylabel("amplitude (arb. unit)", fontsize=fs);
    xticks(fontsize=fs); yticks(fontsize=fs);
    xlim(-0.1,5.1);
    yscale("log"); legend(fontsize=16, loc="center right"); show(); 
