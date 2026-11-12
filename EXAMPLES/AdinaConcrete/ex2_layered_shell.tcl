# AdinaConcrete2D example 2: in-plane loaded wall with layered shells.
# AdinaConcrete2DPlateFiber is used directly as a LayeredShell layer. Units: N, mm, MPa.
wipe
model basic -ndm 3 -ndf 6
set B 1000.0; set H 1000.0; set n 4
for {set j 0} {$j<=$n} {incr j} {for {set i 0} {$i<=$n} {incr i} {
  node [expr {$j*($n+1)+$i+1}] [expr {$i*$B/$n}] [expr {$j*$H/$n}] 0.0}}
for {set i 1} {$i<=$n+1} {incr i} {fix $i 1 1 1 1 1 1}
nDMaterial AdinaConcrete2DPlateFiber 1 30000.0 0.2 2.5 -30.0 -0.002 -6.0 -0.006
uniaxialMaterial Steel01 2 400.0 200000.0 0.01
nDMaterial PlateRebar 3 2 90
nDMaterial PlateRebar 4 2 0
# 150 mm wall: cover/rebar/core/rebar/cover (rebar layers 0.5% each direction, two faces)
section LayeredShell 1 9 1 20 3 0.375 4 0.375 1 26.875 1 26.875 1 26.875 4 0.375 3 0.375 1 46.875
set e 0
for {set j 0} {$j<$n} {incr j} {for {set i 0} {$i<$n} {incr i} {
  set n1 [expr {$j*($n+1)+$i+1}]
  element ShellNLDKGQ [incr e] $n1 [expr {$n1+1}] [expr {$n1+$n+2}] [expr {$n1+$n+1}] 1}}
set top [expr {$n*($n+1)+1}]
for {set i 1} {$i<=$n} {incr i} {equalDOF $top [expr {$top+$i}] 1}
timeSeries Linear 1
pattern Plain 1 1 {load $top 1.0 0 0 0 0 0}
constraints Transformation; numberer RCM; system UmfPack
test NormUnbalance 10.0 100 0
algorithm Newton
integrator DisplacementControl $top 1 0.05
analysis Static
set f [open ex2_curve.txt w]
set ok 0
while {$ok==0 && [nodeDisp $top 1]<3.0} {
  set ok [analyze 1]
  if {$ok!=0} {algorithm KrylovNewton -iterate initial -increment initial; set ok [analyze 1]; algorithm Newton}
  if {$ok==0} {puts $f "[nodeDisp $top 1] [getLoadFactor 1]"}
}
close $f
puts "ex2 finished: top [nodeDisp $top 1] mm, shear [getLoadFactor 1] N, status $ok"
puts "crack state of element 1, point 1, layer 1: [eleResponse 1 material 1 fiber 1 crackState]"
