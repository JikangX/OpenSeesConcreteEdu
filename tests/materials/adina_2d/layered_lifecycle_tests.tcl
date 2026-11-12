# RC membrane coupons with FREE transverse displacements, plastic Steel02
# layers, actual FileDatastore restart, and a mesh comparison. Units: m, MPa.
# Usage: OpenSees.exe layered_lifecycle_tests.tcl native|plain outputDirectory
# "plain" must use the separately built, unchanged forumat.for baseline.
set mode [lindex $argv 0]
set out [lindex $argv 1]
set materialOptions [lindex $argv 2]
if {$mode ni {native plain} || $out eq ""} {error "native|plain and output directory required"}
file mkdir $out
set checks 0
proc require {ok message} {global checks;incr checks;if {!$ok} {error $message}}
proc near {x y tolerance message} {
    require [expr {abs($x-$y)<=$tolerance*max(1,abs($x),abs($y))}] "$message: $x versus $y"
}
proc vectorNear {a b tolerance message} {
    require [expr {[llength $a]==[llength $b]}] "$message dimensions"
    foreach x $a y $b {near $x $y $tolerance $message}
}
proc configure {} {
    wipeAnalysis
    constraints Transformation
    numberer RCM
    system FullGeneral
    test NormUnbalance 1e-9 40 0
    algorithm Newton
    integrator LoadControl 1
    analysis Static
}
proc buildCoupon {n wrapped} {
    global mode control elements nodes materialOptions
    wipe
    model BasicBuilder -ndm 3 -ndf 6
    if {$mode eq "native"} {
        nDMaterial AdinaConcrete2D 1 30000 0 3 -30 -.002 -6 -.006 {*}$materialOptions
        if {$wrapped} {nDMaterial PlateFromPlaneStress 4 1 15000;set concrete 4} else {set concrete 1}
    } else {
        nDMaterial PlaneStressUserMaterial 1 40 7 30 3 -6 -.002 -.006 .001 .05
        nDMaterial PlateFromPlaneStress 4 1 15000
        set concrete 4
    }
    uniaxialMaterial Steel02 2 400 200000 .01 18 .925 .15
    nDMaterial PlateRebar 3 2 0
    section LayeredShell 1 5 $concrete .0396 3 .0006 $concrete .0396 3 .0006 $concrete .0396
    set nodes {};set elements {}
    set control [expr {$n+1}]
    for {set j 0} {$j<=$n} {incr j} {
        for {set i 0} {$i<=$n} {incr i} {
            set tag [expr {1+$i+($n+1)*$j}]
            node $tag [expr {double($i)/$n}] [expr {double($j)/$n}] 0
            lappend nodes $tag
            # Plane membrane motion: ux and uy remain free except supports.
            # Only node 1's uy is fixed, eliminating vertical rigid translation.
            fix $tag [expr {$i==0}] [expr {$tag==1}] 1 1 1 1
            if {$i==$n && $j>0} {equalDOF $control $tag 1}
        }
    }
    set tag 0
    for {set j 0} {$j<$n} {incr j} {
        for {set i 0} {$i<$n} {incr i} {
            incr tag
            set a [expr {1+$i+($n+1)*$j}]
            element ShellMITC4 $tag $a [expr {$a+1}] [expr {$a+$n+2}] [expr {$a+$n+1}] 1
            lappend elements $tag
        }
    }
    timeSeries Linear 1
    pattern Plain 1 1 "load $control 1 0 0 0 0 0"
    configure
}
proc snapshot {} {
    global elements nodes
    set answer [list [getTime] [getLoadFactor 1]]
    foreach nd $nodes {lappend answer [nodeDisp $nd]}
    foreach ele $elements {
        lappend answer [eleResponse $ele forces]
        for {set gp 1} {$gp<=4} {incr gp} {
            lappend answer [eleResponse $ele section $gp deformation] [eleResponse $ele section $gp force]
            for {set layer 1} {$layer<=5} {incr layer} {
                lappend answer [eleResponse $ele section $gp fiber $layer strain] [eleResponse $ele section $gp fiber $layer stress]
            }
        }
    }
    return $answer
}
proc compareSnapshot {a b message} {
    vectorNear [lrange $a 0 1] [lrange $b 0 1] 2e-9 $message
    set index 2
    foreach x [lrange $a 2 end] y [lrange $b 2 end] {
        vectorNear $x $y 2e-9 "$message snapshot-component-$index";incr index
    }
}
proc step {increment {depth 0}} {
    global control elements nodes
    integrator DisplacementControl $control 1 $increment
    set accepted 0
    foreach method {{Newton} {NewtonLineSearch -type Bisection} {KrylovNewton -iterate initial -increment initial}} {
        algorithm {*}$method
        if {[analyze 1]==0} {set accepted 1;break}
    }
    if {!$accepted && $depth<12} {
        # A rejected trial is rolled back by StaticAnalysis. Resolve the SAME
        # requested displacement path in two halves; do not skip an increment
        # or relax force equilibrium to pass a brittle transition.
        step [expr {$increment/2.0}] [expr {$depth+1}]
        step [expr {$increment/2.0}] [expr {$depth+1}]
        return
    }
    require $accepted "RC coupon failed after bounded solver and step refinement at displacement [nodeDisp $control 1]"
    set before [snapshot]
    reactions
    compareSnapshot [snapshot] $before "reaction query changed committed response"
    set residualSquared 0;set support 0;set rightEdgeResidual 0
    foreach nd $nodes {
        set x [nodeCoord $nd 1]
        if {$x==0} {set support [expr {$support+[nodeReaction $nd 1]}]}
        if {$x==1} {
            # equalDOF condenses these individual nodal residuals into ONE
            # retained equation. Audit T^T R: constraint transfer forces at
            # slave nodes need not vanish individually.
            set rightEdgeResidual [expr {$rightEdgeResidual+[nodeReaction $nd 1]}]
        } elseif {$x>0} {
            # A displacement-controlled equation still balances its solved
            # external load; it therefore belongs in the free residual audit.
            set residualSquared [expr {$residualSquared+pow([nodeReaction $nd 1],2)}]
        }
        if {$nd!=1} {set residualSquared [expr {$residualSquared+pow([nodeReaction $nd 2],2)}]}
        near [nodeDisp $nd 1] [expr {$x*[nodeDisp $control 1]}] 2e-9 "uniform membrane strain"
        near [nodeDisp $nd 2] 0 2e-9 "free transverse displacement at nu=0"
    }
    set residualSquared [expr {$residualSquared+$rightEdgeResidual*$rightEdgeResidual}]
    near [expr {sqrt($residualSquared)}] 0 1e-8 "post-commit equilibrium residual"
    near [expr {$support+[getLoadFactor 1]}] 0 1e-8 "global equilibrium"
}
proc continuation {} {
    foreach segment {{100 -.00001} {50 -.00003} {120 -.000005} {150 .00002}} {
        lassign $segment count increment
        for {set k 0} {$k<$count} {incr k} {step $increment}
    }
    return [snapshot]
}
set peaks {}
set wrappers [expr {$mode eq "native" ? {0 1} : {1}}]
foreach wrapped $wrappers {
    foreach n {1 2} {
        puts "CHECK $mode wrapped=$wrapped mesh=$n"
        buildCoupon $n $wrapped
        for {set k 1} {$k<=250} {incr k} {
            step .00001
            if {$k==5} {
                # Independent elastic continuum + steel resultant for unit width.
                near [getLoadFactor 1] [expr {(30000*.1188+200000*.0012)*.00005}] 2e-8 "analytic elastic stiffness"
            }
        }
        set peak [getLoadFactor 1];lappend peaks $peak
        if {$mode eq "native"} {
            # At epsilon=.0025 the original ADINA law releases all tensile
            # concrete stress. Only the two Steel02 layers carry tension.
            near $peak [expr {.0012*[lindex [eleResponse 1 section 1 fiber 2 stress] 0]}] 2e-9 "cracked coupon steel resultant"
        }
        set checkpoint [snapshot]
        set databasePath [file join $out "${mode}_${wrapped}_mesh${n}"]
        database File $databasePath
        save 9
        set reference [continuation]
        wipe
        model BasicBuilder -ndm 3 -ndf 6
        database File $databasePath
        restore 9
        set restored [snapshot]
        foreach kind {checkpoint restored} {
            set f [open "${databasePath}_${kind}.txt" w]
            foreach row [set $kind] {puts $f $row}
            close $f
        }
        compareSnapshot $restored $checkpoint "immediate empty-domain checkpoint restore"
        configure
        compareSnapshot [continuation] $reference "plastic reversal / closure / reopening after restart"
    }
}
foreach peak $peaks {near $peak [lindex $peaks 0] 2e-8 "mesh / direct-versus-wrapped peak invariance"}
puts "PASS $checks $mode RC membrane free-DOF / plastic steel / mesh / FileDatastore lifecycle checks"
wipe
