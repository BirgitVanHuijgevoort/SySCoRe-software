function [PolytopeLarge,PolytopeSmall] = IncreaseDecreasePolytope(Polytope, epsil)
%IncreaseDecreasePolytope Map polytope to a smaller and larger polytope
%based on the bounded error \|y-\hat y\| < eps. 

aux = cell(size(Polytope.A,2),1); 
aux(:) = deal({[-epsil,epsil]});
vertices = combvec(aux{:});
Pol_eps = Polyhedron(vertices');
PolytopeLarge = Polytope+Pol_eps;
PolytopeLarge.minVRep;
PolytopeLarge.computeHRep;

PolytopeSmall = Polytope-Pol_eps;
PolytopeSmall.minVRep;
PolytopeSmall.computeHRep;

end
