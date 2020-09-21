***,H2O
basis=VQZ(f/p)
R=0.95 ANG,THETA=104 DEGREE
geometry={O;H1,O,R;H2,O,R,H1,THETA}
hf              !do closed-shell SCF
