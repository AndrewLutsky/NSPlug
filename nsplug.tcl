##
##NSPlug 1.0
##
##
## Author: Andrew Lutsky
## Use NonStandard Solvents for Solvation
## Inputs required pdbfile of single solvent molecule, psf file of single solvent molecule, density, and size

## Tell Tcl that we're a package and any dependencies we may have
package provide nsplug 1.0

namespace eval ::nsplug:: {
  namespace export nsplug
}

#
# Create the window and initialize data structures
#
#
#NSPlug Tk Window
proc nsplug_tk {} {
  nsplug::nsp
  return $nsplug::w
}


proc nsplug::nsp {} {
  variable w;
  variable webpage "https://www.andrewlutsky.com/";
  variable pdbfile;
  variable psffile;
  variable topfile;
  variable output_path;

  #GENERATION PARAMS
  variable boxsize  0; #Size of Box to Generate
  variable noSim    0; #Don't Run a simulation
  variable shortSim 0; #Run Short Simulation
  variable medSim   0; #Run Med Simulation
  variable longSim  0; #Run Long Simulation
  variable numMols  0; #Number of Molecules to create solution
  
  if { [winfo exists .nsplug] } {
  	wm deiconify $w
  	return
  }
  #set Window options 
  set w [toplevel ".nsplug"]
  wm title $w "nsplug"
  wm resizable $w 1 1
  wm minsize $w 200 200

  #set menubar
  frame $w.menubar -relief raised -bd 2
  pack $w.menubar -fill x
  
  menubutton $w.menubar.help -text "Help" -underline 1 -menu $w.menubar.help.menu -padx 10
  pack $w.menubar.help -side right -fill y
  menu $w.menubar.help.menu -tearoff no
  $w.menubar.help.menu add command -label "About" -command nsplug::about

  #File Loader
  labelframe $w.loadfile -text "Load Files To GENERATE NSS" -relief raised -bd 1
  pack $w.loadfile -side top -fill x -anchor nw -pady 2

  frame $w.loadfile.buttons 
  pack $w.loadfile.buttons -side left -anchor nw -fill y -expand 1 -padx 10 -pady 10

  button $w.loadfile.buttons.loadpdb -text "Load PDB" -command nsplug::loadpdb
  pack $w.loadfile.buttons.loadpdb -side left -anchor nw -fill y -expand 1 -padx 10	
 
  #button $w.loadfile.buttons.loadpsf -text "Load PSF" -command nsplug::loadpsf
  #pack $w.loadfile.buttons.loadpsf -side left -anchor n -fill y -expand 1 -padx 10
  
  button $w.loadfile.buttons.loadtop -text "Load Topology" -command nsplug::loadtop
  pack $w.loadfile.buttons.loadtop -side right -anchor ne -fill y -expand 1 -padx 10

  #Solvent Generator Parameter
  labelframe $w.gen -text "NSS Generation Parameters" -relief raised -bd 1
  pack $w.gen -side top -fill x -anchor nw -pady 2


  #BoxSize Parameter
  frame $w.gen.boxsize
  pack $w.gen.boxsize -side left -fill x -anchor nw -pady 2 -padx 2
  label $w.gen.boxsize.label -text "Box Size of NSS"
  entry $w.gen.boxsize.entry -bg White -bd 2 -width 5 -textvariable nsplug::boxsize
  pack $w.gen.boxsize.label $w.gen.boxsize.entry -side left -anchor nw -pady 2
 
  #MaxCount Parameter 
  #frame $w.gen.density
  #pack $w.gen.density -side left -fill x -anchor nw -pady 2 -padx 2
  #label $w.gen.density.label -text "Number of Molecules to Pack"
  #entry $w.gen.density.entry -bg White -bd 2 -width 5 -textvariable nsplug::numMols
  #pack $w.gen.density.label $w.gen.density.entry -side left -anchor nw -pady 2



  #Output Path
  frame $w.out
  pack $w.out -side top -anchor nw -fill x -pady 2 -padx 2
  label $w.out.lab -text "Output Path"
  entry $w.out.val -bg White -bd 2 -width 20 -textvariable nsplug::output_path
  pack $w.out.lab $w.out.val -side left -anchor nw -pady 2

  #Create the NSS
  frame $w.go
  pack $w.go -side bottom -fill x -anchor sw -pady 2 -padx 2
  button $w.go.createnss -text "Create NSS" -command nsplug::createnssLat
  pack $w.go.createnss -side bottom -fill x -anchor sw -pady 2 -padx 2

}


proc nsplug::about {} {
  puts " Made by Andrew Lutsky : www.andrewlutsky.com "


}

proc nsplug::loadpdb {} {
  variable pdbfile
  set pdbfile [tk_getOpenFile -title "pdb" -filetypes [list {"PDB files" {.pdb}}] ]
  if {![file readable $pdbfile]} {
    puts "File Not Readable"
  }
  puts $pdbfile
}


#Deprecated
#proc nsplug::loadpsf {} {
#  variable psffile 
#  set psffile [tk_getOpenFile -title "psf" -filetypes [list {"PSF files" {.psf}}]]
#
#  if {![file readable $psffile]} {
#    puts "File Not Readable"
#  }
#  puts $psffile
#}


proc nsplug::loadtop {} {
  variable topfile 
  set topfile [tk_getOpenFile -multiple true -title "top" -filetypes [list {"TOP files" {.rtf .str}}]]

  if {![file readable $topfile]} {
    puts "File Not Readable"
  }
  puts $topfile

}

#Housekeeping Functions to find the max radius and center the solvent molecule
proc nsplug::maxsize { sel } {
  set mm [ measure minmax $sel]
  set xdiff [ expr { [lindex $mm 1 0] - [lindex $mm 0 0] } ]
  set ydiff [ expr { [lindex $mm 1 1] - [lindex $mm 0 1] } ]
  set zdiff [ expr { [lindex $mm 1 2] - [lindex $mm 0 2] } ]
  
  set lis "$xdiff $ydiff $zdiff"
  set dis [ expr { ($xdiff ** 2 + $ydiff ** 2 + $zdiff ** 2) ** 0.5 } ]
  set rad [ expr { $dis / 2. } ]
  return $rad
}

proc nsplug::center {sel} {
  set cen [measure center $sel]
  $sel moveby [vecscale -1 $cen]
}


proc nsplug::genpsf {numcopies} {
  variable topfile
  package require psfgen 2.0
  resetpsf


  psfgen_logfile "NSSPSF.log"
  foreach a [split $topfile " "] {
  	topology $a
  }
  #I suppose if you want to solvate with proteins?
  pdbalias atom ILE CD1 CD ; # formerly "alias atom ..."
  pdbalias residue HIS HSE ; #aliases residue HIS to HSE
  
  segment NSS {   	
    for {set i 1} {$i <= $numcopies} {incr i} {
      pdb outp$i.pdb
    }
  }
  for {set i 1} {$i <= $numcopies} {incr i} {   
  	coordpdb outp$i.pdb NSS 
  }

  guesscoord
  writepsf output.psf
  writepdb output.pdb

  psfgen_logfile close
}

proc nsplug::createnssLat {} {
  draw material "Ghost"
  variable output_path
  variable pdbfile
  variable boxsize
  #dummy pdb file to set radius
  mol new $pdbfile waitfor all
  set sel [atomselect top all]
  nsplug::center $sel
  set radius [nsplug::maxsize $sel]
  mol delete top
  
  puts "Working for you"

  cd $output_path
  file mkdir output
  cd output

  set sum 0;
  for {set k 0} { $k < $boxsize} {incr k} {
    for {set j 0} {$j < $boxsize} {incr j} {
      for {set i 0} {$i < $boxsize} {incr i} {
	set x [expr { (2 * $i + (($j + $k) % 2)) * $radius } ]
        set y [expr { sqrt(3) * $radius * ($j + (1. / 3.) * ($k % 2)) } ]
	set z [expr { ( (2. / 3.) * sqrt(6) * $radius * $k) } ]
        
	#loads in new molecule and centers it at 0 0 0 
	mol new $pdbfile waitfor all
	set sel [atomselect top all]
        nsplug::center $sel	
	
	#rotates and moves molecule
	set randomrotx [ expr { rand() * 360 } ]
	set randomroty [ expr { rand() * 360 } ]
	set randomrotz [ expr { rand() * 360 } ]
	$sel move [transaxis x $randomrotx]
	$sel move [transaxis y $randomroty]
	$sel move [transaxis z $randomrotz]
	$sel moveby "$x $y $z"
	
	incr sum;
	$sel set resid $sum
	$sel writepdb outp$sum.pdb
        
	mol delete top
	#draw sphere "$x $y $z" radius "$radius" resolution 40
        unset x
        unset y
        unset z
        unset randomrotx
        unset randomroty
        unset randomrotz


      }
    }
  }
  
  nsplug::genpsf $sum
  unset sum
  mol new output.psf
  mol addfile output.pdb  
}



