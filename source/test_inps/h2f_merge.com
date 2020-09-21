***,example for merge
print,orbitals,basis
rh2=1.4
rhf=300.
basis=vdz
symmetry,x,y                            !use C2v symmetry
geometry={F}

text,F
{rhf;wf,9,1,1;occ,3,1,1;orbital,2130.2} !rhf for f-atom

text,H2
symmetry,x,y                            !use C2v symmetry
geometry={
         H1,
         H2,H1,rh2}

{hf;orbital,2100.2}                     !scf for h2
{multi;occ,2;orbital,2101.2}            !mcscf for h2

text,FH2
geometry={F;                            !linear geometry for F+H2
         H1,F,rhf
         H2,H1,rh2,F,180}

{merge
orbital,2130.2                          !rhf orbitals for F-atom
move,1.1,2.1,1.1                        !move orbitals 1.1, 2.1
move,3.1,0.4,4.1;                       !move all remaining, starting at 4.1
orbital,2100.2                          !hf orbitals for H2
move,1.1,0.4                            !move these to free positions
save,2131.2}                            !save merged orbitals

{rhf;occ,4,1,1;start,2131.2             !rhf for F+H2
orbital,2132.2}

{merge
orbital,2130.2                          !rhf orbitals for F-atom
move,1.1,2.1,1.1                        !move orbitals 1.1, 2.1
move,3.1,3.1,4.1;                       !move orbital 3.1 to 4.1
move,4.1,0.4,6.1                        !move all remaining, starting at 6.1
orbital,2101.2                          !mcscf orbitals for H2
move,1.1,0.4                            !move these to free positions
save,2141.2}                            !save merged orbitals

{multi;occ,5,1,1;start,2141.2}          !casscf for F+H2 using valence space
