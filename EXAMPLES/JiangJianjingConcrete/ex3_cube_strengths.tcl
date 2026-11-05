# JiangJianjingConcrete example 3: stress-controlled unit cube (one stdBrick),
# uniaxial, biaxial and confined compression strengths (lecture 7, p. 77:
# MSC.Marc gives 30 / 36 / 44 / 66.6 MPa for fc = 30 MPa).
# Stage 1 applies the lateral confinement, stage 2 increases the active
# stresses by 0.05 MPa per step until the solution fails; the last converged
# value is the strength. Units: N, mm, MPa.
proc cube {} {
    wipe
    model basic -ndm 3 -ndf 3
    nDMaterial JiangJianjingConcrete 1 30.0 30000.0 0.18
    set k 0
    foreach z {0.0 1.0} {foreach {x y} {0 0 1 0 1 1 0 1} {node [incr k] $x $y $z}}
    fix 1 1 1 1; fix 2 0 1 1; fix 3 0 0 1; fix 4 1 0 1; fix 5 1 1 0; fix 6 0 1 0; fix 8 1 0 0
    element stdBrick 1 1 2 3 4 5 6 7 8 1
}
proc faceLoads {sx sy sz} {
    foreach n {2 3 6 7} {load $n [expr {$sx/4.0}] 0 0}
    foreach n {3 4 7 8} {load $n 0 [expr {$sy/4.0}] 0}
    foreach n {5 6 7 8} {load $n 0 0 [expr {$sz/4.0}]}
}
proc solver {} {
    constraints Plain; numberer Plain; system FullGeneral
    test NormUnbalance 1e-9 400 0
    algorithm ModifiedNewton -initial  ;# or: -tangent secant + Newton
}
proc strength {lateral active} {
    cube
    solver
    if {$lateral != 0.0} {
        timeSeries Linear 1
        pattern Plain 1 1 "faceLoads [expr {-$lateral}] [expr {-$lateral}] 0.0"
        integrator LoadControl 0.05; analysis Static
        for {set i 0} {$i < 20} {incr i} {analyze 1}
        loadConst -time 0.0
        wipeAnalysis
        solver
    }
    timeSeries Linear 2
    pattern Plain 2 2 "faceLoads {*}$active"
    integrator LoadControl 0.05; analysis Static
    set last 0.0
    for {set i 0} {$i < 3000} {incr i} {
        if {[analyze 1] != 0} break
        set last [getLoadFactor 2]
    }
    return $last
}
puts "ex3 uniaxial            [strength 0.0 {0.0 0.0 -1.0}] MPa (Marc 30)"
puts "ex3 biaxial s2=s3       [strength 0.0 {0.0 -1.0 -1.0}] MPa (Marc 36)"
puts "ex3 confinement 2 MPa   [strength 2.0 {0.0 0.0 -1.0}] MPa (Marc 44)"
puts "ex3 confinement 6 MPa   [strength 6.0 {0.0 0.0 -1.0}] MPa (Marc 66.6)"
