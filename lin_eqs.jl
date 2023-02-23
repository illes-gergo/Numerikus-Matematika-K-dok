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
    A2[m, n] .-= A2[m, k] ./ A2[k, k] .* A1[k, n]
    b2[m] .-= A2[m, k] ./ A2[k, k] * b2[k]
end

x2 = zeros(size(A1, 2), 1)

for m in size(A1,1):(-1):1
    local temp = 0
end