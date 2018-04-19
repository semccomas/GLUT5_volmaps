# we seem to also need pbc wrapping here as well

mol new /data2/GLUT5/POPC_POPE.90_10/production/fc0_bb/PC_PE_fc0.1.gro
animate read xtc /data2/GLUT5/POPC_POPE.90_10/production/fc0_bb/PC_PE_fc0.1-10.skip50.xtc waitfor all 0

pbc wrap -centersel "name BB" -compound res -all

set f [molinfo 0 get numframes] 
set compare [atomselect 0 "name BB"]
set reference [atomselect 0 "name BB" frame 0]
set ref_minmax [measure minmax [atomselect 0 "name BB" frame 0]]
set compare_system [atomselect 0 "all"]
for {set i 0} {$i < $f} {incr i} {
    #rename for correct frame
    $compare frame $i
    # compute transformation
    set trans_mat [measure fit $compare $reference]
    $compare_system frame $i
    #do alignment
    $compare_system move $trans_mat    
    volmap occupancy [atomselect 0 "name BB SC1 SC2 SC3 SC4 and within 6 of resname POPE" frame $i] -res 0.5 -mol 0 -o $i.POPE.dx -minmax $ref_minmax
}

