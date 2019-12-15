% Funcion para puntos de desplazamiento
function z = centers(theta,m,r)
    z = r*exp(i*2*theta*m);
end
