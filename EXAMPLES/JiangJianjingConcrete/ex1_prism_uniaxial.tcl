# JiangJianjingConcrete example 1: plain-concrete prism (2 stdBrick), uniaxial
# compression by prescribed top displacement. Units: N, mm, MPa.
# Shows: the equivalent uniaxial curve (Chinese-code parabola, peak fc at
# eps0 = 2 fc / E0), the stress hold after beta = 1, and crushing.
wipe
model basic -ndm 3 -ndf 3
set a 100.0; set h 200.0
set k 0
foreach z [list 0.0 [expr {$h/2}] $h] {foreach {x y} [list 0 0 $a 0 $a $a 0 $a] {node [incr k] $x $y $z}}
fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1
nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18
element stdBrick 1 1 2 3 4 5 6 7 8 1
element stdBrick 2 5 6 7 8 9 10 11 12 1
equalDOF 9 10 3; equalDOF 9 11 3; equalDOF 9 12 3
timeSeries Linear 1
pattern Plain 1 1 {load 9 0 0 -1.0}
constraints Transformation; numberer Plain; system FullGeneral
test NormUnbalance 1e-3 400 0
algorithm ModifiedNewton -initial
integrator DisplacementControl 9 3 -0.01
analysis Static
set f [open ex1_curve.txt w]
puts $f "axial_strain axial_stress_MPa beta crushed"
for {set i 1} {$i <= 120} {incr i} {
  if {[analyze 1] != 0} {puts "stopped at step $i"; break}
  set st [eleResponse 1 material 1 state]
  puts $f "[expr {-[nodeDisp 9 3]/$h}] [expr {-[getLoadFactor 1]/($a*$a)}] [lindex $st 1] [lindex $st 0]"
}
close $f
puts "ex1 finished: strain [expr {-[nodeDisp 9 3]/$h}], stress [expr {-[getLoadFactor 1]/($a*$a)}] MPa"
