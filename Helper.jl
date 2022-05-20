using LinearAlgebra
using Plots
plotly()


n = 11
x_domain = [-1. 1.;
		    -1. 1.]

h = 2 / (n -1)

a(x) = [-x[2],x[1]]
# a(x) = [x[1],x[2]]

u_0(x) =  norm(x,2) < 1 ? exp(-1/(1-norm(x,2)^2)) : 0
u_0(x) = 1 


function analytical_solution(t,x)
	Δt = 0.01
	path = Matrix{Float64}(undef,2,1) 
	path[1] = x[1] 
	path[2] = x[2] 

	while t > Δt
		x = x .- Δt .* a(x)
		path = hcat(path,x)
		t = t - Δt
	end

	x = x .- t .* a(x)	
	path = hcat(path,x)
	return (u_0(x),x, path) 
end

function plot_path(sol)
	plot(sol[3][1,:], sol[3][2,:], xlim=[-1,1], ylim=[-1,1], label="")
end

sol = analytical_solution(3,[-1,-1])
plot_path(sol)

# Parameter
begin
	n = 5
	x_domain = [-1. 1.;
			    -1. 1.]
	
	h = 2 / (n -1) 

	a(x) = [-x[2],x[1]]
	u_0(x) =  norm(x,2) < 1 ? exp(-1/(1-norm(x,2)^2)) : 0
	u_0(x1_, x2_) =  sqrt(x1_^2 + x2_^2)^2 < 1 ? exp(-1/(1-sqrt(x1_^2 + x2_^2)^2)) : 0

	x1 = LinRange(x_domain[1,1],x_domain[1,2],n)
	x2 = LinRange(x_domain[2,1], x_domain[2,2],n)
end

t = 0
σ = 1/h

begin
	mutable struct Blob
		x::Vector{Number}
		m::Number
	end


	function evolve_blobs(t::Number, blobs::Vector{Blob})
		Δt = 0.01

		while t > Δt
			for i in 1:size(blobs)[1]
					blobs[i].x = blobs[i].x + Δt * a(blobs[i].x)
					t = t -Δt
					@info t
			end
		end

		for i in 1:size(blobs)[1]
				blobs[i].x = blobs[i].x + Δt * a(blobs[i].x)
		end
	end


	ζ(x) = 1/(2π) *	(6-6*norm(x,2)^2 + norm(x,2)^4) * exp(-norm(x,2)^2)

	σζ(x) = σ^-2 * ζ(x./σ)

	x = LinRange(-2,2,100)
	plot(x, σζ.(x))



	function reconstruct_gird_values(blobs::Vector{Blob})
		u = zeros(n,n)
		
		for i in 1:n
			for j in 1:n
				for k in 1:size(blobs)[1]
					u[i,j] += blobs[k].m * σζ( blobs[k].x - [x1[i],x2[j]])
				end
			end
		end
		
		return u
	end




	function plot_blobs(t_end)

		blobs = [Blob([i,j] + h/2 .*[1,1], h^2 * u_0([i,j])) for i in x1[1:end] for j in x2[1:end]]
		@info blobs[1].x

		evolve_blobs(t_end, blobs)


		ζ(x) = 1/(2π) *	(6-6*norm(x,2)^2 + norm(x,2)^4) * exp(-norm(x,2)^2)
		
		σζ(x) = σ^-2 * ζ(x./σ)


		plot_ = heatmap(reconstruct_gird_values(blobs))

		for i in 1:size(blobs)[1]
				scatter!(plot_, (blobs[i].x[1] , blobs[i].x[2]), label="", c="red")
				# scatter!(plot_, (blobs[i].x[1] * 1/h + n/2, blobs[i].x[2] * 1/h + n/2), label="", c="red", markeralpha=0.3, markersize=20)
		end

		return plot_
	end
end
plot_blobs(t)

