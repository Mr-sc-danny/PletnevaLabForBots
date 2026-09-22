###############################################################
# save_hemo_frames_by_folder.tcl                                #
# DESCRIPTION:                                                  #
#   For each pdb/dcd system below, loads the structure +        #
#   trajectory, then writes ONE pdb per frame containing only   #
#   the "resname HEMO" atoms, into that system's own folder     #
#   under ruffling_MD (matching the existing folder layout).    #
#                                                                 #
#   Each written frame is then immediately repaired in place    #
#   (see fix_hemo_pdb below) so it conforms to the strict        #
#   80-column PDB spec that PorphyStruct/ChemSharp expects.      #
#   VMD's writepdb only writes out to column 78 (it never adds  #
#   the optional charge field, cols 79-80), and the "HEMO"       #
#   residue name is 4 letters where the standard code is the     #
#   3-letter "HEM", so the residue-name/chain-ID fields end up   #
#   glued together with no separator. Both issues make          #
#   PorphyStruct throw "Specified argument was out of range of  #
#   valid values" when it tries to load these files.             #
#                                                                 #
# USAGE (from the Tk console):                                   #
#         source save_hemo_frames_by_folder.tcl                  #
#                                                                 #
# NOTE: Edit the $systems list below if the folder-to-system     #
#       mapping (cryt vs NMR) does not match what you intended.  #
###############################################################

# ---------------------------------------------------------------
# Rebuilds a single ATOM/HETATM line into a proper, full 80-column
# PDB record:
#   - pads the line out to 80 columns (adds the missing charge field)
#   - collapses any 4+ letter "HEM*" residue name down to the
#     standard 3-letter "HEM", restoring the blank separator
#     column before the chain ID
#   - fills in the element symbol (cols 77-78) when VMD left it
#     blank (this happens for hydrogens), inferring it from the
#     first alphabetic character of the atom name
# ---------------------------------------------------------------
proc fix_atom_line {line} {
    if {[string length $line] < 80} {
        append line [string repeat " " [expr {80 - [string length $line]}]]
    }

    set rec      [string range $line 0 5]
    set serial   [string range $line 6 10]
    set atomname [string range $line 12 15]
    set altloc   [string index $line 16]
    set resname  [string trim [string range $line 17 19]]
    if {[string match "HEM*" $resname]} { set resname "HEM" }
    set chain    [string index $line 21]
    set resseq   [string range $line 22 25]
    set icode    [string index $line 26]
    set x        [string range $line 30 37]
    set y        [string range $line 38 45]
    set z        [string range $line 46 53]
    set occ      [string range $line 54 59]
    set temp     [string range $line 60 65]
    set elem     [string trim [string range $line 76 77]]

    if {$elem eq ""} {
        set elem ""
        foreach ch [split $atomname ""] {
            if {[string is alpha -strict $ch]} {
                set elem [string toupper $ch]
                break
            }
        }
    }

    return [format "%-6s%5s %-4s%1s%-3s %1s%4s%1s   %8s%8s%8s%6s%6s          %2s  " \
        $rec $serial $atomname $altloc $resname $chain $resseq $icode \
        $x $y $z $occ $temp $elem]
}

# ---------------------------------------------------------------
# Reads a just-written pdb file and rewrites it in place with
# every ATOM/HETATM line passed through fix_atom_line. Other
# record types (CRYST1, END, etc.) are left untouched.
# ---------------------------------------------------------------
proc fix_hemo_pdb {path} {
    set fin [open $path r]
    set content [read $fin]
    close $fin

    set fout [open $path w]
    foreach line [split $content "\n"] {
        set line [string trimright $line "\r"]
        if {[string length $line] == 0} { continue }
        set rec [string range $line 0 5]
        if {$rec eq "ATOM  " || $rec eq "HETATM"} {
            puts $fout [fix_atom_line $line]
        } else {
            puts $fout $line
        }
    }
    close $fout
}

# ---- Source directory (existing pdb/dcd files) ----
set dist_meas_dir "E:/Desktop/research/MD/slit/dist_meas"

# ---- Destination base directory (existing folder layout) ----
set ruffling_dir "E:/Desktop/research/MD/slit/ruffling_MD"

# ---- Systems to process: {output_folder pdb_file dcd_file} ----
set systems {
    {1YNR "1ynr.pdb"  "1ynr.dcd"}
}

foreach entry $systems {
    set folder  [lindex $entry 0]
    set pdbfile [lindex $entry 1]
    set dcdfile [lindex $entry 2]

    set pdbpath  "$dist_meas_dir/$pdbfile"
    set dcdpath  "$dist_meas_dir/$dcdfile"
    set outdir   "$ruffling_dir/$folder"

    puts "\n==============================================="
    puts "System: $folder"
    puts "  pdb: $pdbpath"
    puts "  dcd: $dcdpath"
    puts "  out: $outdir"
    puts "==============================================="

    # Make sure the destination folder exists
    file mkdir $outdir

    # Load structure + trajectory into a fresh molecule
    set molid [mol new $pdbpath waitfor all]
    mol addfile $dcdpath type dcd waitfor all molid $molid

    set n [molinfo $molid get numframes]
    set hemoSel [atomselect $molid "resname HEMO"]

    set start_frame 1000
    set end_frame   3000
    for {set i $start_frame} {$i <= $end_frame && $i < $n} {incr i} {
        $hemoSel frame $i
        set outfile "$outdir/HEMO_${i}.pdb"
        $hemoSel writepdb $outfile
        fix_hemo_pdb $outfile
        puts "\t progress ($folder): $i/$n"
    }

    puts "\t progress ($folder): $end_frame/$end_frame -- Done."
    puts "\t wrote frames $start_frame-$end_frame (clamped to $n total) to $outdir"

    $hemoSel delete
    mol delete $molid
}

puts "\nAll systems processed."