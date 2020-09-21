***,h2o population analysis
geometry={o;h1,o,r;h2,o,r,h1,theta}   !Z-matrix geometry input
r=1 ang                               !bond length
theta=104                             !bond angle
basis=6-311g**
hf                                    !do scf calculation
pop;                                  !Mulliken population analysis using mcscf density
individual                            !give occupations of individual basis functions
