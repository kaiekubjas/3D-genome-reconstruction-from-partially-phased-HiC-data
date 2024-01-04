using HomotopyContinuation, Random, Combinatorics

# Fix random number generator seed for reproducibility
Random.seed!(123)

# Variables and parameters
@var x[1:3,1:12], y[1:3,1:12], C[1:12,1:12]

# The original rational functions from system (3.7) in the paper
normsq = v->sum([v[i]^2 for i=1:length(v)])

full_system = [C[i,j]-1/normsq(x[:,i]-x[:,j])-1/normsq(x[:,i]-y[:,j])-1/normsq(y[:,i]-x[:,j])-1/normsq(y[:,i]-y[:,j]) for i=1:12 for j=i+1:12]

# Up to the action of O(3,C)\rtimes C^3, we can assume that the following coordinates are zero
reduction = [x[1,1],x[2,1],x[3,1],x[1,2],x[2,2],x[1,3]]

# Form a parametric square system of rational functions
param = [C[i,j] for i=1:12 for j=i+1:12]
F = System(vcat(full_system,reduction),parameters=param)

# Create a random start solution for the monodromy
x0 = randn(3,12)
y0 = randn(3,12)
x0[1,1],x0[2,1],x0[3,1],x0[1,2],x0[2,2],x0[1,3] = zeros(6)
sol0=vcat(x0[:],y0[:])

# Compute the corresponding parameter values
param0 = [1/normsq(x0[:,i]-x0[:,j])+1/normsq(x0[:,i]-y0[:,j])+1/normsq(y0[:,i]-x0[:,j])+1/normsq(y0[:,i]-y0[:,j]) for i=1:12 for j=i+1:12]

# Check that we really have a (very good approximation of a) solution to the system for this choice of parameters
# evaluate(F(sol0),param=>param0)

# Function that for a given solution computes the orbit under the aciton of (Z/2)^10 on the solutions
# Note that we broke the x_i<->y_i symmetry for i=1,2,3 when we partly fixed the coordinates of x_1,x_2,x_3
# However, we also have an symmetry (x,y)<->(-x,-y)
function symmetry(v)
    xpart = reshape(v[1:36],3,12)
    ypart = reshape(v[37:72],3,12)
    
    all_solutions = []
    
    for subset in powerset(4:12)
        xpartnew = deepcopy(xpart)
        ypartnew = deepcopy(ypart)
        for I in subset
                xpartnew[:,I] = deepcopy(ypart[:,I])
                ypartnew[:,I] = deepcopy(xpart[:,I])
        end
        append!(all_solutions,[vcat(xpartnew[:],ypartnew[:]),
                -vcat(xpartnew[:],ypartnew[:])])
    end
    return all_solutions
end

# Check that all points in the orbit of sol0 are (approximate) solutions
# [norm(evaluate(F(sol),param=>param0)) for sol in symmetry(sol0)]

# Solve using monodromy (stop after 1001 solutions are found)
MS = monodromy_solve(F, sol0, param0, group_action=symmetry, target_solutions_count=1001)

# Certify
C = certify(F,solutions(MS),target_parameters=param0)

# Save the certificates as a text file
save("ambiguous_ID_certificates.txt",C)
