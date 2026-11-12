# AdinaConcrete2D example 1: plane-stress RC cantilever strip (quad + truss rebar),
# displacement-controlled push. Units: N, mm, MPa.
wipe
model basic -ndm 2 -ndf 2
set L 1000.0; set H 200.0; set t 150.0; set nx 10; set ny 4
for {set j 0} {$j<=$ny} {incr j} {
  for {set i 0} {$i<=$nx} {incr i} {
    node [expr {$j*($nx+1)+$i+1}] [expr {$i*$L/$nx}] [expr {$j*$H/$ny}]
  }
}
for {set j 0} {$j<=$ny} {incr j} {fix [expr {$j*($nx+1)+1}] 1 1}
# Concrete: E nu ft fc epsc fu epsu (compression negative). Defaults = ADINA program.
nDMaterial AdinaConcrete2D 1 30000.0 0.2 2.5 -30.0 -0.002 -6.0 -0.006
set e 0
for {set j 0} {$j<$ny} {incr j} {
  for {set i 0} {$i<$nx} {incr i} {
    set n1 [expr {$j*($nx+1)+$i+1}]
    element quad [incr e] $n1 [expr {$n1+1}] [expr {$n1+$nx+2}] [expr {$n1+$nx+1}] $t PlaneStress 1
  }
}
# Bottom and top bars (area 400 mm2) along the first and last node rows.
uniaxialMaterial Steel01 2 400.0 200000.0 0.01
for {set i 0} {$i<$nx} {incr i} {
  element truss [incr e] [expr {$i+1}] [expr {$i+2}] 400.0 2
  set top [expr {$ny*($nx+1)+$i+1}]
  element truss [incr e] $top [expr {$top+1}] 400.0 2
}
set tip [expr {$ny*($nx+1)+$nx+1}]
timeSeries Linear 1
pattern Plain 1 1 {load $tip 0.0 -1.0}
constraints Transformation; numberer RCM; system UmfPack
test NormUnbalance 1.0 100 0
algorithm Newton
integrator DisplacementControl $tip 2 -0.05
analysis Static
set f [open ex1_curve.txt w]
set ok 0
while {$ok==0 && [nodeDisp $tip 2]>-8.0} {
  set ok [analyze 1]
  if {$ok!=0} {algorithm KrylovNewton -iterate initial -increment initial; set ok [analyze 1]; algorithm Newton}
  if {$ok==0} {puts $f "[expr {-[nodeDisp $tip 2]}] [getLoadFactor 1]"}
}
close $f
puts "ex1 finished: tip [nodeDisp $tip 2] mm, load [getLoadFactor 1] N, status $ok"
