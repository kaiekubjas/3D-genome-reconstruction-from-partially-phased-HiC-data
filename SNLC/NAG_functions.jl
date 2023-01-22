using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, Random, IterTools, Plots, Statistics

# Function that takes as input a list of lists [L1,...,LN] where Li is a list of points in R^n, and outputs a 
# tensor T\in R^(L1\times ...\times LN) where T(v1,...,vN)=\sum_{i=1}^N ||v_i-(v1+...+vN)/N||.
function dist_tensor(listlist)
    N = length(listlist)
    D=zeros([length(L) for L in listlist]...)
    for idxtuple in Iterators.product([1:length(L) for L in listlist]...);
        center = 1/N*sum([listlist[i][idxtuple[i]] for i=1:N])
        D[idxtuple...] = sum([norm(listlist[i][idxtuple[i]]-center) for i=1:N])/N
    end
    return D
end;      

# Function that takes a vector v in R^n as input, and outputs ||v||^2.
normsq = v->sum([v[i]^2 for i=1:length(v)]);

# Functions that create the square polynomial system from the paper.
make_polynomial = (x,y,z,c) -> ( c*normsq(z-x)*normsq(z-y)-normsq(z-x)-normsq(z-y) );
make_system = (x,y,Z,C) -> [make_polynomial(x,y,Z[:,i],C[i]) for i=1:6];

# Function that takes a list of points in C^6, and remove duplicates up to the action (a1,a2,a3,a4,a5,a6) -> (a4,a5,a6,a1,a2,a3).
function one_solution_per_orbit(solution_list)
    reduced_solution_list = Vector{Vector{ComplexF64}}(undef,0)
    for s in solution_list
        sprime = vcat(s[4:6],s[1:3])
        add_to_list = true
        for t in reduced_solution_list
            if norm(t-sprime)<1e-10
                add_to_list = false
                break
            end
        end
        if add_to_list
            append!(reduced_solution_list,[s])
        end
    end
    return reduced_solution_list
end;

# Function that takes a list of points in C^6, and for each point (a1,a2,a3,a4,a5,a6), it adds (a4,a5,a6,a1,a2,a3) to the list.
function union_of_orbits(reduced_solution_list)
    solution_list = Vector{Vector{ComplexF64}}(undef,0)
    for s in reduced_solution_list
        append!(solution_list,[s,vcat(s[4:6],s[1:3])])
    end
    return solution_list
end;

# Estimation of the ambiguous loci using homotopy continuation
function estimate_ambig_htpy(Pred,Zua,num_ua_sets=5;param0=nothing,result0=nothing,number_of_contacts=20,max_attempts=100,imag_part_threshold=0.1)
    
    m = size(Pred,2) # number of ambiguous loci
    n = floor(Int64,m+size(Pred,1)/2) #total nuber of loci
    
    # Create the general system
    @var x[1:3];
    @var y[1:3];
    @var Z[1:3,1:6];
    @var C[1:6];
    F = System([make_polynomial(x,y,Z[:,i],C[i]) for i=1:6],parameters=vcat(Z[:],C))
    
    # If we haven't already provided a solution set for a certain set of parameters, 
    # then we just make a random (complex) choice of parameters, and solve with polyhedral start system
    if isnothing(param0) || isnothing(result0)
        # Solve the system for a random choice of complex parameters
        Z0 = randn(ComplexF64,3,6)
        C0 = randn(ComplexF64, 6)
        param0 = vcat(Z0[:],C0)
        result0 = solve(F, target_parameters = param0)
    end
    
    # Initialize array in which we will put the estimations given as output
    Xambig = Array{Float64}(undef,3,m)
    Yambig = Array{Float64}(undef,3,m)

    # Initialize array which we will use to keep track of the number of 
    # attempts needed for each bead, before we find real solutions
    # 1st coordinate: index of the set of unambiguous beads
    # 2nd coordinate: index of the amboguous bead
    attempts = zeros(Int64,num_ua_sets,m)

    # Initialize array which we will use to store the choices of unambiguous beads
    ua_choices = zeros(Int64,m,num_ua_sets,6)
    
    # Initialize array which we will use to store the minimal distances
    minimal_distances = Array{Float64}(undef,m)

    start = time()

    for j in 1:m
        
        print("Locus: ")
        print(j)
        print(" of ")
        print(m)
        print("\n")
        

        # Initialize array for storing the lists of real solutions
        realsols = Array{Vector{Vector{Float64}}}(undef,num_ua_sets)
        
        # Pick out the relevant column of partially unambiguous contacts
        P_column = Pred[:,j]
        
        # Pick out indices of the the highest contacts counts
        P_column = replace!(P_column,NaN=>-Inf)
        sorted_P_column_indices = sort(1:2*(n-m),rev=true,by=i->P_column[i])
        ua_indices_highest_P = sorted_P_column_indices[1:number_of_contacts]

        for k in 1:num_ua_sets

            realsols[k] = []
            
            # Keep making attempts until we find real solutions (maximally 300 attempts)
            while length(realsols[k])==0 && attempts[k,j]<max_attempts
                
                #print("\n")
                
                # Increase the attempt counter by one
                attempts[k,j] = attempts[k,j]+1
                
                # Make a choice of six unambiguous beads (redo it if it's not unique)
                ua_choices[j,k,:] = sort(ua_indices_highest_P[randperm(number_of_contacts)[1:6]])
                while ua_choices[j,k,:] in [ua_choices[j,i,:] for i=1:k-1] || length(filter(y->y<=n-m,ua_choices[j,k,:]))<2 || length(filter(y->y>n-m,ua_choices[j,k,:]))<2
                    ua_choices[j,k,:] = sort(randperm(2*(n-m))[1:6])
                end
                
                sol = solve(
                     F,
                     one_solution_per_orbit(solutions(result0));
                     start_parameters =  param0,
                     target_parameters = vcat(Zua[:,ua_choices[j,k,:]][:],Pred[ua_choices[j,k,:],j])
                 )
                
                 # Pick out the real solutions. If no solution is found within max_attempts
                 # we instead use an average of the neighbouring unambiguous beads. 
                 if attempts[k,j]<max_attempts
                     realsols[k] = union_of_orbits(real_solutions(sol,real_tol=imag_part_threshold))
                 else
                    realsols[k] = [vcat(mean(Zua[:,filter(y->(y<=n-m),ua_choices[j,k,:])],dims=2),
                        mean(Zua[:,filter(y->(y>n-m),ua_choices[j,k,:])],dims=2))[:]]
                    print("Estimated by averaging.\n")
                 end
            end
            
        end

        # Compute a distance tensor, and find which solutions from the tuple of runs are the most similar
        D = dist_tensor([realsols[i] for i=1:num_ua_sets])
        minimal_distances[j] = findmin(D)[1]
        minidx = findmin(D)[2]
        

        # Compute their mean
        meansolution = 1/num_ua_sets*sum([realsols[i][minidx[i]] for i=1:num_ua_sets])
                
        # Take the resulting coordinates as our estimation
        Xambig[:,j] = meansolution[1:3]
        Yambig[:,j] = meansolution[4:6]

    end

    T = time()-start;
    
    return Xambig, Yambig, T, attempts, ua_choices, minimal_distances
	
end;
