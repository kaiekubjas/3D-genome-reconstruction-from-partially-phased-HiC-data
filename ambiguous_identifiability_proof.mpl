restart;
RandomTools[MersenneTwister][SetState](123);


# Number of loci
n := 12:

# Form the map psi
normsq := (x,y,z) -> x^2+y^2+z^2:
for i from 1 to n do
	for j from i+1 to n do
		C[i,j] := 1/normsq((x[i,1],x[i,2],x[i,3])-(x[j,1],x[j,2],x[j,3])) 
		+ 1/normsq((x[i,1],x[i,2],x[i,3])-(y[j,1],y[j,2],y[j,3])) 
		+ 1/normsq((y[i,1],y[i,2],y[i,3])-(x[j,1],x[j,2],x[j,3])) 
		+ 1/normsq((y[i,1],y[i,2],y[i,3])-(y[j,1],y[j,2],y[j,3]));
	end do:
end do:

psi := [seq(seq(C[i,j],j=(i+1)..n),i=1..n)]:

# Jacobian
Jac := Student[VectorCalculus][Jacobian](psi,[seq(seq(x[i,j],j=1..3),i=1..n),seq(seq(y[i,j],j=1..3),i=1..n)]):

# Pick a random point with integer coordinates
max_integer := 100:
for i from 1 to n do
	for j from 1 to 3 do
		xspec[i,j]:=rand(1..max_integer)():
		yspec[i,j]:=rand(1..max_integer)():
	end do
end do;

# Check that xspec and yspec belong to the affine open subset where psi is defined
for i from 1 to n do
	for j from i+1 to n do
		if normsq((xspec[i,1],xspec[i,2],xspec[i,3])-(xspec[j,1],xspec[j,2],xspec[j,3]))*normsq((xspec[i,1],xspec[i,2],xspec[i,3])-(yspec[j,1],yspec[j,2],yspec[j,3]))*normsq((yspec[i,1],yspec[i,2],yspec[i,3])-(xspec[j,1],xspec[j,2],xspec[j,3]))*normsq((yspec[i,1],yspec[i,2],yspec[i,3])-(yspec[j,1],yspec[j,2],yspec[j,3])) = 0 then
			print("Not generic enough!");
		end if
	end do
end do;

# Evaluate the Jacobian
Jac_evaluated := subs([seq(seq(x[i,j]=xspec[i,j],j=1..3),i=1..n),seq(seq(y[i,j]=yspec[i,j],j=1..3),i=1..n)],Jac):

# Compute the rank
LinearAlgebra[Rank](Jac_evaluated);


