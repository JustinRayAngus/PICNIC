sample inputs for neutrals driven by specular reflection off
of a piston moving at a consant velocity. The neutrals are atomic
deuterium and the initial density is spatially uniform at 1e20/m^3

piston_collisionless - collisionless neutrals initilized with
                       zero energy in the piston direction

piston_collisional - same initialization as piston_collisionless,
                     but VHS collisions are turned on

piston_vp1e2 - neutral atomic deuterium gas initialized with 300 K 
               temperature (VT ~= 1 km/s) being compressed by piston 
               moving at Vp = 100 m/s (Vp<<VT ==> adiabiatic compression)

piston_vp1e3 - neutral atomic deuterium gas initialized with 300 K 
               temperature (VT ~= 1 km/s) being compressed by piston 
               moving at Vp = 1000 m/s (Vp~=VT ==> weak shock)

piston_vp1e4 - neutral atomic deuterium gas initialized with 300 K 
               temperature (VT ~= 1 km/s) being compressed by piston 
               moving at Vp = 10000 m/s (Vp>>VT ==> strong shock)
