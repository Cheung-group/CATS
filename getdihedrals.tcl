
set chan [open dihedralangles.out w+]
set frames [molinfo top get numframes]

for {set i 0} {$i < $frames} {incr i} {
set sel [atomselect top "name CA" frame $i]
set dihedrals [$sel get {phi psi}]
puts $chan "$dihedrals"

}

close $chan