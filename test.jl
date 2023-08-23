using Plots

Lx = 1
Ly = 1
k = 4*pi/Lx
l = 4*pi/Ly

uᵢ(x, y) = -(1/(k^2+l^2))*sin(k*x)*(l*cos(l*y)-32/Ly^2*y*sin(l*y)+16/Ly*sin(l*y))*exp(-16/Ly^2*y^2+16/Ly*y-4)
vᵢ(x, y) = k/(k^2+l^2)*cos(k*x)*sin(l*y)*exp(-16/Ly^2*y^2+16/Ly*y-4)

#default(size=(600,600), fc=:heat)
#=
z = []
for x in LinRange(0,1,11)
    for y in LinRange(0,1,11)
        print(uᵢ(x,y))
    end
end
=#

X = LinRange(0,1,11)
Y = LinRange(0,1,11)

#surface(X,Y,(X, Y)->uᵢ(X,Y), linealpha = 0.3)
surface(X,Y,(X, Y)->vᵢ(X,Y), linealpha = 0.3)