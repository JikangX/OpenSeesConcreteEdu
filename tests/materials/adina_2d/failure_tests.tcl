# End-to-end rejection and recovery: an excessive material increment must
# abort the global step, preserve the committed shell state, and permit a
# smaller retry. A separate spring supplies a nonempty system of equations.
proc require {condition message} {if {!$condition} {error $message}}
proc same {a b message} {
    require [expr {[llength $a]==[llength $b]}] "$message size"
    foreach x $a y $b {require [expr {abs($x-$y)<1e-12*max(1,abs($x),abs($y))}] $message}
}
wipe
set shellType [lindex $argv 0]
if {$shellType eq ""} {set shellType ShellMITC4}
if {$shellType ni {ShellMITC4 ShellNLDKGQ}} {error "unsupported shell type"}
model BasicBuilder -ndm 3 -ndf 6
nDMaterial AdinaConcrete2D 1 30000 .2 3 -30 -.002 -6 -.006 -maxStrainIncrement 2e-5
section LayeredShell 1 3 1 .12
node 1 0 0 0;node 2 1 0 0;node 3 1 1 0;node 4 0 1 0
element $shellType 1 1 2 3 4 1
node 101 2 0 0;node 102 2 0 0
fix 101 1 1 1 1 1 1;fix 102 0 1 1 1 1 1
uniaxialMaterial Elastic 2 1000
element zeroLength 2 101 102 -mat 2 -dir 1
timeSeries Path 1 -dt 1 -values 0 1e-8 .5 2e-8
pattern Plain 1 1 {
    foreach nd {1 2 3 4} {
        sp $nd 1 [nodeCoord $nd 1]
        foreach dof {2 3 4 5 6} {sp $nd $dof 0}
    }
}
timeSeries Linear 2
pattern Plain 2 2 {load 102 1 0 0 0 0 0}
constraints Transformation
numberer RCM
system FullGeneral
test NormUnbalance 1e-12 10
algorithm Newton
integrator LoadControl 1
analysis Static
require [expr {[analyze 1]==0}] "seed shell step"
set force [eleResponse 1 forces]
set stress [eleResponse 1 material 1 fiber 1 stress]
set strain [eleResponse 1 material 1 fiber 1 strain]
require [expr {[llength $strain]==5}] "material response available"
require [expr {[analyze 1]!=0}] "constitutive rejection was swallowed by $shellType"
require [expr {abs([getTime]-1)<1e-12}] "failed step did not revert domain time"
same [eleResponse 1 material 1 fiber 1 stress] $stress "failed step stress rollback"
same [eleResponse 1 material 1 fiber 1 strain] $strain "failed step strain rollback"
same [eleResponse 1 forces] $force "failed step internal force rollback"
# Force/stiffness/reaction queries must leave the constitutive history intact.
reactions
same [eleResponse 1 material 1 fiber 1 stress] $stress "reactions changed committed stress"
integrator LoadControl 1e-6
require [expr {[analyze 1]==0}] "smaller retry after section rejection"
require [expr {[getTime]>1}] "recovery did not advance"
puts "PASS shell constitutive failure propagation / rollback / query purity / recovery"
wipe
