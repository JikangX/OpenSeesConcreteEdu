# AdinaConcrete3D example 3: plain-concrete prism, uniaxial compression by
# prescribed top displacement, lateral faces free. Units: N, mm, MPa.
wipe
model basic -ndm 3 -ndf 3
set a 100.0; set h 200.0
set k 0
foreach z [list 0.0 [expr {$h/2}] $h] {foreach {x y} [list 0 0 $a 0 $a $a 0 $a] {node [incr k] $x $y $z}}
# symmetry supports at the base
fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1
nDMaterial AdinaConcrete3D 1 30000.0 0.2 2.5 -30.0 -0.002 -6.0 -0.006
element stdBrick 1 1 2 3 4 5 6 7 8 1
element stdBrick 2 5 6 7 8 9 10 11 12 1
equalDOF 9 10 3; equalDOF 9 11 3; equalDOF 9 12 3
timeSeries Linear 1
pattern Plain 1 1 {load 9 0 0 -1.0}
constraints Transformation; numberer Plain; system FullGeneral
test NormUnbalance 1.0 50 0
algorithm Newton
integrator DisplacementControl 9 3 -0.01
analysis Static
set f [open ex3_curve.txt w]
for {set i 1} {$i<=40} {incr i} {
  if {[analyze 1]!=0} {puts "stopped at step $i"; break}
  puts $f "[expr {-[nodeDisp 9 3]/$h}] [expr {[getLoadFactor 1]/($a*$a)}]"
}
close $f
puts "ex3 finished: strain [expr {-[nodeDisp 9 3]/$h}], stress [expr {[getLoadFactor 1]/($a*$a)}] MPa"
