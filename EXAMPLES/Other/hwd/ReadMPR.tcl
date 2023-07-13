## -----------------------------------------------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  Date: January 19, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# READ MPR FILE FOR SECTION PROPERTIES SCRIPT

# CSDir 	-> 	The directory containing the mpr files with the cross section information, e.g. "./Dimensions/DeckCS/"
# CSType 	->	The type of cross section (Deck or Cap)
# BentNum 	-> 	The bent number corresponding to the cross section (for deck, the 1st bent or abutment out of the 2) (ib)
# NodeNum 	->  (For deck) The intermediate node number corresponding to the cross section, the 1st node out of the 2 (iNode)
# NR  		->  (For cap beam) Set to 1 if the cap beam is on the NR side, default 0

# execute as: ReadMPR $CSDir $CSType $BentNum {-NodeNum $NodeNum -NR $NR}

proc ReadMPR {CSDir CSType BentNum optArgs} {
    # Set the defaults for NodeNum and NR
    array set options {-NodeNum 0 -NR 0};
    # Read in the NodeNum and NR arguments
    array set options $optArgs;
	set NodeNum $options(-NodeNum);
	set NR $options(-NR);
	
	global in; # Read in unit and constant variables
	
	##### Get the file corresponding to the specified cross section
	if {$CSType == "Deck"} {
		# If it's a deck cross section
		if {$BentNum <= 5 || ($BentNum == 6 && $NodeNum == 0)} {
			set fileName [glob -path $CSDir *Abut1-Bent6*];					# The cross sections from Abut1.0-Bent6.0 are all the same
		} elseif {$BentNum == 8 && $NodeNum <= 3} {
			set fileName [glob -path $CSDir *Bent8-8-9_3*]; 				# The cross sections from Bent8.0-Bent8.3 are all the same
		} elseif {$BentNum <= 7 || $BentNum <= 11} {
			set fileName [glob -path $CSDir *Bent$BentNum*_$NodeNum*]; 		# For Bent6.1-Bent7.4 and Bent8.4-Bent11.4, naming is Bent(ib)-(ib+1)_(iNode)
		} elseif {$NodeNum <= 4} { 											# NE side takes nodes 0-4
			set fileName [glob -path $CSDir *Bent$BentNum*_$NodeNum*_NE*]; 	# For Bent12.0-Bent14.4 (NE side), naming is Bent(ib)-(ib+1)_(iNode)_NE
		} elseif {$BentNum >= 14 || ($BentNum >= 13 && $NodeNum >=8)} {		# NR side takes nodes 5-9.
			set fileName [glob -path $CSDir *13-14_3-Abut15_NR*]; 			# Bent13.8-14.9 (NR side) are all the same
		} else { 															
			set NodeNumNR [expr $NodeNum-5];
			set fileName [glob -path $CSDir *Bent$BentNum*_$NodeNumNR*_NR*];# For Bent 12.5-13.7 (NR side), naming is Bent(ib)-(ib+1)_(iNode-5)_NR
		}
	} else {
		# If it's a cap beam cross section
		if {$NR} {  # Bents 13 and 14 have a different section for the NR side
			puts "NRside"
			set fileName [glob -path $CSDir *B$BentNum-NR-f.mpr];
			# set fileName [glob -path $CSDir *B$BentNum-NR.mpr];			
		} else {
		set fileName [glob -path $CSDir *B$BentNum-f.mpr];
		# set fileName [glob -path $CSDir *B$BentNum.mpr];		
		}
	}
	# puts "$fileName";  # Check the file name
	
	##### Read the Area, MOIz, MOIy, and J from the MPR file
	set fpCS [open $fileName r];							# Open the file, get the filepath
	set CSID [regexp -inline {[^\/]+$} $fileName]; 			# Read the identifying description of the filename (e.g. "Bent8-9_4.mpr")
	# puts "CSID = $CSID" 	# Check the file description
	set Nline 0;
	foreach line [split [read $fpCS] \n] {
		incr Nline; 										# Start reading each line; count how many lines have been read
		if {$Nline == 4} {
			set wordloop 0;
			foreach word [split $line { }] { # Split by spaces; read each word and count how many words have been read
				if {$word !=""} {
					incr wordloop;
					if {$wordloop == 2} {
						set A [expr $word*$in**2];			# Get the Area, the 2nd word of the 4th line
						# puts "file=$CSID, Nline=$Nline, wordloop=$wordloop A=$A"
					}
				}
			}
		} elseif {$Nline == 16} {
			set wordloop 0;			
			foreach word [split $line { }] { # Split by spaces; read each word and count how many words have been read
				if {$word !=""} {
					incr wordloop;
					if {$wordloop == 2} {
						set Iy [expr $word*$in**4];			# Get the MOI about the horizontal (local y) axis, the 2nd word of the 16th line
						# puts "file=$CSID, Nline=$Nline, wordloop=$wordloop, Iy=$Iy"
					}
				}
			}
		} elseif {$Nline == 17} {
			set wordloop 0;
			foreach word [split $line { }] { # Split by spaces; read each word and count how many words have been read
				if {$word !=""} {
					incr wordloop;
					if {$wordloop == 2} {
						set Iz [expr $word*$in**4];			# Get the MOI about the vertical (local z) axis, the 2nd word of the 17th line
						# puts "file=$CSID, Nline=$Nline, wordloop=$wordloop, Iz=$Iz"
					}
				}
			}
		}
	}
	set J [expr $Iy+$Iz];
	# puts "file=$CSID, J=$J";
	close $fpCS;
	# puts "ib=$BentNum, iNode=$NodeNum, A=$A, Iz=$Iz, Iy=$Iy, J=$J"
	return "$A $Iy $Iz $J"; # Return a list of the section properties
}