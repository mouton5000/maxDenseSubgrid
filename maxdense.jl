using Distributions

using JuMP

#using GLPK

#using CoinOptServices

using CPLEX

# Definiing parameter

#######################################################################################################################################




# number of rows


m= int(ARGS[1])

# number of columns
n= int(ARGS[2])

# density distributions
p=0.1

#we define the grid
G=zeros(Int64,m,n)
for i=1:m
	for j=1:n
		if (rand()<p)
			G[i,j]=1
		end
	end
end
println("matrice G= ",G)

# We define the first column addition matrix (C_{n-1})
C=eye(Int64,n)
C[n,n]=0
C[n,n-1]=1
#println(C)
println("")

# We define the first lign addition matrix (C_{m-1})
L=eye(Int64,m)
L[m,m]=0
L[m-1,m]=1
#println(L)
println("")

T=m+n-1
#####################################################################################################################################################################################




#Definition du model 
mod = Model(solver=CplexSolver())





# M[i,j,t] matrix 
@defVar(mod,  M[i=1:m,j=1:n,t=1:T] >= 0)

# x variables
@defVar(mod, 0 <= x[t=1:T] <= 1, Int )

# produit M[i,j,t]x[t]
@defVar(mod,  mx[i=1:m,j=1:n,t=1:T]  >= 0)

# produit M[i,j,T]*M[k,j,T]
@defVar(mod,  mm[i=1:m,j=1:n,k=1:m,l=1:n]  >= 0)


# completement stupide ... la prochaine fois ajouter lignes et colonnes de zeros
@setObjective(mod, Max, 0.5*(mm[1,1,1,2]+mm[1,1,2,1]+mm[1,1,2,2]+mm[m,1,m,2]+mm[m,1,m-1,1]+mm[m,1,m-1,2]+ mm[1,n,1,n-1]+mm[1,n,2,n] + mm[1,n,2,n-1]+mm[m,n,m,n-1]+mm[m,n,m-1,n]+mm[m,n,m-1,n-1]+sum{mm[i,1,i-1,1]+mm[i,1,i+1,1]+mm[i,1,i,2]+mm[i,1,i-1,2]+mm[i,1,i+1,2], i=2:(m-1)} + sum{mm[1,i,1,i-1]+mm[1,i,1,i+1]+mm[1,i,2,i]+mm[1,i,2,i-1]+mm[1,i,2,i+1], i=2:(n-1)}+ sum{mm[i,n,i-1,n]+mm[i,n,i+1,n]+mm[i,n,i,n-1]+mm[i,n,i-1,n-1]+mm[i,n,i+1,n-1], i=2:(m-1)} + sum{mm[m,i,m,i-1]+mm[m,i,m,i+1]+mm[m,i,m-1,i]+mm[m,i,m-1,i-1]+mm[m,i,m-1,i+1], i=2:(n-1)}+sum{mm[i,j,i-1,j-1]+mm[i,j,i,j-1]+mm[i,j,i+1,j-1]+mm[i,j,i+1,j]+mm[i,j,i+1,j+1]+mm[i,j,i,j+1]+mm[i,j,i-1,j+1]+mm[i,j,i-1,j], i=2:(m-1), j=2:(n-1)}))


# M_1=G
for i=1:m
	for j=1:n
		@addConstraint(mod, M[i,j,1]==G[i,j])
	end
end

for t=2:n
	Ct=C-eye(Int64,n)
	for i=1:m
		for j=1:n
			@addConstraint(mod, M[i,j,t]==M[i,j,t-1]+sum{Ct[k,j]*mx[i,k,t-1], k=1:n})
		end
	end
	# mise a de C
	if(t<n)
		C[n-t+1,n-t+1]=0
		C[n-t+1,n-t]=1
	end
#	println(C)
end

for t=n+1:T
	Lt=L-eye(Int64,m)
	for i=1:m
		for j=1:n
			@addConstraint(mod, M[i,j,t]==M[i,j,t-1]+sum{Lt[i,k]*mx[k,j,t-1], k=1:m})
		end
	end
	# mise a de C
	if(t<T)
		ind=t-n
		L[m-ind,m-ind]=0
		L[m-ind-1,m-ind]=1
	end
#	println(L)
end

#linearisation
for i=1:m
	for j=1:n
		@addConstraint(mod, M[i,j,T] <= 1 )
		for t=1:T
			@addConstraint(mod, mx[i,j,t] <= M[i,j,t] )
			@addConstraint(mod, mx[i,j,t] <= x[t] )
			@addConstraint(mod, mx[i,j,t] >= M[i,j,t]+x[t]-1 )
		end
		for k=1:m
			for l=1:n
				@addConstraint(mod, mm[i,j,k,l] <= M[i,j,T] )
				@addConstraint(mod, mm[i,j,k,l] <= M[k,l,T] )
				@addConstraint(mod, mm[i,j,k,l] >= M[i,j,T]+M[k,l,T]-1 )
			end
		end
	end
end


# On commence par resoudre la relaxation P^0 ou il n'ya aucune contrainte induite par le second niveau
#print(mod)
tic()
status = solve(mod)
println("valeur_optimale ",getObjectiveValue(mod))
toc()
mat=getValue(M)
for t=1:T
	println("iteration ",t)
	println("matrice= ",round(Int64,mat[:,:,t]))
	println("")
end


   

