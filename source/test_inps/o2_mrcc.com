***,O2 tests
memory,8,m
gthresh,energy=1.d-8

geometry={o1;o2,o1,r1}
r1=2.2
set,state=1,symmetry=4,spin=2  ! Triplet sigma- state
basis=vdz

rhf
uccsd(t)
method(1)='UCCSD(T) MOLPRO'
e(1)=energy

rccsd(t)
method(2)='RCCSD(T) MOLPRO'
e(2)=energy

mrcc,method=ccsdt,dir=mrccdir
method(3)='CCSDT MRCC'
e(3)=energy

mrcc,method=ccsdtq,restart=1,dir=mrccdir
method(4)='CCSDT MRCC'
e(4)=energy

table,method,e
