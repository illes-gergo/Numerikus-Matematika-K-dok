using LinearAlgebra;

# 2. Solution ha M == N
A = [2; 4;; 3; 5];
b = [-2; 7];

x = A \ b
println(x)

x = A^(-1) * b
println(x)

x = inv(A) * b
println(x)

# Mátrix szinguláris -> det(A) = 0
try # -> próbáld meg
    local A = [3; 6;; 5; 10]
    local b = [2; 5]
    local x = A \ b
catch # -> hiba van
    println("Nem sikerült!")
    det([3; 6;; 5; 10])
end

# M < N -> alulhatározott egyenletrendszer

A = [1; 2;; 2; 3;; 3; 4]
b = [4; 5]

alpha = inv(A * transpose(A)) * b

xplus = transpose(A) * alpha
println(xplus)
x = pinv(A) * b
println(x)

# M > N -> túlhatározott egyenlet

# Gauss ellimináció

A0 = [1 -1 3; 2 -1 2; 2 -1 1]
b0 = [0, -4, -7]
x0 = A0 \ b0

# Előre

A1 = copy(A0)
b1 = copy(b0)

for k in 1:(size(A1, 1)-1), m in (k+1):size(A1, 1), n in (k+1):size(A1, 1)
    A1[m, n] -= A1[m, k] / A1[k, k] * A1[k, n]
    b1[m] -= A1[m, k] / A1[k, k] * b1[k]
end


x1 = zeros(size(A1, 2), 1)

for m in size(A1, 1):(-1):1
    local temp = 0
    for n in (m+1):size(A1, 1)
        temp += A1[m, n] * x1[n]
    end
    x1[m] = (b1[m] - temp) / A1[m, m]
end

# Gauss ellimináció mátrixműveletekkel

A2 = copy(A0)
b2 = copy(b0)

for k in 1:(size(A2, 1)-1)
    m = (k+1):size(A2, 1)
    n = (k+1):size(A2, 1)
    A2[m, n] .-= A2[m, k] / A2[k, k] * A2[k, n]'
    b2[m] .-= A2[m, k] / A2[k, k] * b2[k]
end

x2 = zeros(size(A2, 2), 1)

for m in size(A2, 1):(-1):1
    n = (m+1):size(A2, 1)
    x2[m] = (b2[m] - sum(A2[m, n] .* x2[n])) / A2[m, m]
end

# PLU Felbontás

A00 = [1 2 5; 0.2 1.6 7.4; 0.5 4 8.5]

if size(A00)[1] != size(A00)[2]
    error("A mátrix nem négyzetes!")
end

A = copy(A00)

N = size(A, 1)
P = collect(1.0I(N))
U = zeros(size(A))
L = collect(1.0I(N))
for k in 1:N
    while A[k, k] == 0
        A[[end, k], :] .= A[[k, end], :]
        P[[end, k], :] .= P[[k, end], :]
    end

    for n in k:N
        U[k, n] = A[k, n]
    end

    for m in (k+1):N
        L[m, k] = A[m, k] / A[k, k]
        A[m, k] = L[m, k]
        for n in (k+1):N
            A[m, n] = A[m, n] - A[m, k] * A[k, n]
        end
    end
end

println("L*U - P*A00")
display(L * U - P * A00)

# PLU Felbontás mátrixműveletekkel

if size(A00)[1] != size(A00)[2]
    error("A mátrix nem négyzetes!")
end

A = copy(A00)

N = size(A, 1)

P = collect(1.0I(N))
U = zeros(size(A))
L = collect(1.0I(N))

for k in 1:N
    while A[k, k] == 0
        A[[end, k], :] .= A[[k, end], :]
        P[[end, k], :] .= P[[k, end], :]
    end
    U[k, k:N] .= A[k, k:N]
    L[(k+1):N, k] .= A[(k+1):N, k] / A[k, k]
    A[(k+1):N, k] .= L[(k+1):N, k]
    A[(k+1):N, (k+1):N] .= A[(k+1):N, (k+1):N] .- A[(k+1):N, k] * A[k, (k+1):N]'
end

println("Mátrixokkal L*U - P*A00")
display(L * U - P * A00)