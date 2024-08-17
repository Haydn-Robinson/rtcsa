# ----------------------------------------------------------------------------------------------------------
# Real Time Control Simulation and Analysis
# Author: Haydn Robinson
# Version: 4.0
#
# A library developed for use in simulating and analysing industrial aerospace propulsion control systems given input data.
#
# ----- IMPORTANT INFO ------
# Many functions require the sampling time (same as iteration rate/task period) of the data. This must be defined at the global level as sampTime.
# E.g: set sampTime 0.025
# The functions will call on this global variable instead of requiring the sample time as an input.
#
# If a function requires the previous value of an input or output, it must be initialised in the first pass logic to a sensible value.
# E.g prevInput could be set to the current value of the input for the first pass).

# Contents
# - Validation Utilities
#	- Detect Missing Signals (X)
# 	- Time Step Validation (X)
#	- Channel Change Detection (X)
# 	- Common code for TCL analysis (X)
#	- Enable Output (X)

# - Functions
#	- Load Signal (X)
#	- Maximum value 
#	- Minimum Value
#	- Mean Value (Average)
#	- Moving Average
#	- Boolean Logic Functions
#		- XOR
#	- Protected Divide
#	- Differentiation
#	- Integration
#	- Latch
#	- First Order Lag
#	- Greater Than with Hysteresis
#	- Less Than with Hysteresis
#	- Bracket Lookup Function
#	- 1D Interpolation
#	- 2D Interpolation
#	- 1D Graphic Data Lookup
#	- 2D Graphic Data Lookup
#	- Fault Integrator
#	- Delay On
#	- Delay Off
#	- True/False Detect
#	- False/True Detect
#	- Divide By Zero Protection


namespace eval rtcsa {

# ----------------------------------------------------------------------------------------------------------
# Data Rate Control

# Use this utility to downsample data from base data rate to iteration rate of system being replicated.

proc IterationRateControl {count} {
	global sampTime; global dataRate
	set multiplier [expr {round($sampTime/$dataRate)}];
	if {[expr {$count >= $multiplier}]} {
		set isInvocationPoint 1; set count 0
	} else {
			set isInvocationPoint 0
	}
	incr count 1
	return [list $count $isIterationPoint]
}

# ----------------------------------------------------------------------------------------------------------
# Maximum value function

# This function takes any number of input arguments and returns the maximum. The online EMS script analysis site does not recognise the min/max functions,
# this library uses this as an alternative.

proc MyMax {args} {
	set output -inf
	foreach arg $args {
		if {$arg > $output} {
			set output $arg
		}
	}
	return $output
}
	
# ----------------------------------------------------------------------------------------------------------
# Minimum value function

# This function takes any number of input arguments and returns the minimum. The online EMS script analysis site does not recognise the min/max functions,
# this library uses this as an alternative.

proc MyMin {args} {
	set output inf
	foreach arg $args {
		if {$arg < $output} {
			set output $arg
		}
	}
	return $output
}

# ----------------------------------------------------------------------------------------------------------
# Mean value function

# This function takes any number of input arguments and calculates their mean value.

proc Mean {args} {
	set sum 0
	set length [llength $args]
	foreach arg $args {
		set sum [expr {$sum + $arg}]
		}
	return [expr {$sum/double($length)}] 
}

# ----------------------------------------------------------------------------------------------------------
# Statistical Analysis function

# This function takes a list as an input, and calculates its mean and population standard deviation. The stdType argument chooses between the 
# population standard deviation (default - leave blank or enter '0') and sample standard deviation (enter '1').

proc Statistic {values {stdType 0}} {
	set sum 0; set stdSum 0
	set length [llength $values]
	foreach val $values {
		set sum [expr {$sum + $val}]
		}
	set mean [expr {$sum/double($length)}]
	foreach val $values {
		set stdSum [expr {$stdSum + pow(($val - $mean),2)}]
	}
	if {stdType == 0} {
		set norm $length
	} elseif {stdType == 1} {
		Set norm [expr {$length – 1}]
	}
	set std [expr {sqrt((1/double($norm))*$stdSum)}]
	return [list $mean $std]
}


# ----------------------------------------------------------------------------------------------------------
# Moving Average function

# This function performs a simple, un-weighted moving average on the input data. The windowSize argument defines the size of the window
# (in samples) over which the function will average the data. The signal argument is the current value of the signal being averaged, whilst
# buffer is a list of length windowSize that contains the data samples being used to compute the moving average.

proc MovingAverage {windowSize signal buffer} {
	lappend $buffer $signal
	if {[llength $buffer] > $windowSize} {
		set buffer [lreplace $buffer 0 0]
	}
	set avgVal [Mean $buffer]
	return [list $avgVal $buffer]
}


# ----------------------------------------------------------------------------------------------------------

# Exclusive Or (XOR) // Modulo 2 addition

# 0 + 0 = 0
# 0 + 1 = 1
# 1 + 0 = 1			Where + is Modulo-2 addition
# 1 + 1 = 0

proc XOR {input1 input2} {
	if {[expr {($input1 && !$input2) || ($input2 && !$input1)}]} {
		set output 1
	} else {
		set output 0
	}
	return output
}
# ----------------------------------------------------------------------------------------------------------
# Protected Divide

# Function to perform division with protection from a zero denominator and integer rounding errors. Function also prevents Inf output.
	
proc SafeDivide {numerator denominator} {
	set protectedDenom [expr {double($denominator)}]
	if {$protectedDenom >= 0} {
		set protectedDenom [MyMax $protectedDenom [expr {1.17755e-38}]]
	} else {
		set protectedDenom [MyMin $protectedDenom [expr {-1.17755e-38}]]
	}
	set output [expr {$numerator/$protectedDenom}]
	return $output
}

# ----------------------------------------------------------------------------------------------------------
# Differentiation

# Performs differentiation using Euler backwards differentiation approximation
# Udot(k) = [U(k) - U(k-1)]/T			(T - sampling time)

proc Differentiate {input prevInput} {
	global sampTime
	set output [SafeDivide [expr {$input - $prevInput}] $sampTime]
	return $output
}

# ----------------------------------------------------------------------------------------------------------
# Integration

# Performs integration using Euler backwards differentiation approximation:
# Y(k) = T*U(k) + Y(k-1)				(T - sampling time)
# Output is limited to between integration limits minLim and maxLim and is resettable

proc Integrate {input minLim maxLim prevOutput reset} {
	global sampTime
	set output [expr {$sampTime*$input + $prevOutput}]
	set output [MyMax $minLim $output]; set output [MyMin $maxLim $output]
	if {$reset == 1} {
		set output $input
	}
	return $output
}

# ----------------------------------------------------------------------------------------------------------
# Latch

# If clear is true, the output is false. Else if trigger is true, the output is set true. Otherwise, the output will hold its previous value.
 
proc Latch {trigger clear prevOutput} {
	if {$clear == 1} {
		return 0
	} elseif {$trigger == 1} {
		return 1
	} else {
		return $prevOutput
	}
}
# ----------------------------------------------------------------------------------------------------------
# First Order discrete Lag (Low pass filter)

# Performs a first order lag using Euler backwards differentiation approximation. 
# Y(k) = aU(k) + (1-a)Y(k-1)		where a = T/(T + TLag)

# Y(k) = (T*U(k) + TLag*Y(k-1))/(T+TLag)
                
# T - sampling time
# TLag - lag time constant

proc Lag {lagTimeConst input prevOutput} {
	global sampTime
	set output [SafeDivide [expr {$sampTime*$input + $lagTimeConst*$prevOutput}] [expr {$sampTime + $lagTimeConst}]]
	return $output
}

# ----------------------------------------------------------------------------------------------------------
#Greater than with Hysteresis function.

# Greater than function with hysteresis

proc GTHyst {input setpoint band prevOutput} {
	if {$input > [expr {$setpoint + $band}]} {
		set GTHystOut 1
	} elseif {$input < [expr {$setpoint - $band}]} {
		set GTHystOut 0
	} else {
      set GTHystOut $prevOutput
	}
	return $GTHystOut
}

# ----------------------------------------------------------------------------------------------------------
#Less than with Hysteresis function.

# Less than function with hysteresis

proc LTHyst {input setpoint band prevOutput} {
	if {$input < [expr {$setpoint - $band}]} {
		set LTHystOut 1
	} elseif {$input > [expr {$setpoint + $band}]} {
		set LTHystOut 0
	} else {
      set LTHystOut $prevOutput
   }
	return $LTHystOut
}
# ----------------------------------------------------------------------------------------------------------
# Bracket Lookup Function

# Takes inputs of xVals and x, where xVals is a list of lookup points for the input x. The function returns the upper
# and lower index of the interval in which the value x is found. For use in graphic data lookup functions.

proc BracketLookup {xVals x} {
	set idx1 0; set idx2 [expr {[llength $xVals] - 1}]
	while {[expr {($idx2 - $idx1) > 1}]} {
		set idxMid [expr {($idx2 + $idx1)/2}]
		if {$x > [lindex $xVals $idxMid]} {
			set idx1 $idxMid
		} else {
			set idx2 $idxMid
		}
	}
	return [list $idx1 $idx2]
}

# ----------------------------------------------------------------------------------------------------------
# Linear Interpolation:

# Performs a linear (1D) interpolation of the input. Output is limited to values between fx1 and fx2 inclusive.
# Function determines if x values are increasing or decreasing, and then limits input x to within the range [x1,x2].

# fx2|                o
#    |                               f(x1)     f(x)   f(x2)
# fx |         x                       |--------|------|
#    |                                 x2       x      x2
# fx1|   o
#    |_________________                fx = fx1 + [(fx2 - fx1)(x - x1)]/(x2 - x1)                        
#       x1     x     x2

proc Interpolate1D {x1 x2 fx1 fx2 x} {
	set boundaryUsed 0
	set xDiff [expr {$x2 - $x1}]
	if {$xDiff > 0} {
		if {$x > $x2} {
			set fx $fx2; set boundaryUsed 1
		} elseif {$x < $x1} {
			set fx $fx1; set boundaryUsed 1
		}
	} elseif {$xDiff < 0} {
		if {$x < $x2} {
			set fx $fx2; set boundaryUsed 1
		} elseif {$x > $x1} {
			set fx $fx1; set boundaryUsed 1
		}
	} else {
		set fx $fx2; set boundaryUsed 1
	}
	if {!$boundaryUsed} {
		set xRatio [expr {double(($x - $x1))/$xDiff}]
		set fx [expr {$fx1 + (($fx2 - $fx1)*$xRatio)}]
	}
	return $fx
}

# ----------------------------------------------------------------------------------------------------------
# Bilinear Interpolation:

# Performs a Bilinear (2D) interpolation of the input. Output is limited to the co-domain bordered by fx1y1, fx1y2, fx2y1 and fx2y2 inclusive.
# Function first determines if y values are increasing or decreasing, and then limits input y to within the range [y1,y2]. If the input is in
# range, it performs an interpolation to determine fx1y and fx2y, and then interpolates again using Interpolate1D to find fxy.

proc Interpolate2D {x1 y1 x2 y2 fx1y1 fx1y2 fx2y1 fx2y2 x y} {
	set boundaryUsed 0
	set yDiff [expr {$y2 - $y1}]
	if {$yDiff > 0} {
		if {$y > $y2} {
			set fx1y $fx1y2; set fx2y $fx2y2; set boundaryUsed 1
		} elseif {$y < $y1} {
			set fx1y $fx1y1; set fx2y $fx2y1; set boundaryUsed 1
		}
	} elseif {$yDiff < 0} {
		if {$y < $y2} {
			set fx1y $fx1y2; set fx2y $fx2y2; set boundaryUsed 1
		} elseif {$y > $y1} {
			set fx1y $fx1y1; set fx2y $fx2y1; set boundaryUsed 1
		}
	} else {
		set fx1y $fx1y2; set fx2y $fx2y2; set boundaryUsed 1
	}
	if {!$boundaryUsed} {
		set yRatio [expr {double(($y - $y1))/$yDiff}]
		set fx1y [expr {$fx1y1 + $yRatio*($fx1y2 - $fx1y1)}]
        set fx2y [expr {$fx2y1 + $yRatio*($fx2y2 - $fx2y1)}]
	}
	set fxy [Interpolate1D $x1 $x2 $fx1y $fx2y $x]
	return $fxy
}

# ----------------------------------------------------------------------------------------------------------
# 1D lookup

# Takes inputs of xVals and fxVals, lists of equal length that define the lookup points of the input (x) and the output (fxVals).
# Function performs a one dimensional graphical lookup using the input x, determining the output by interpolation.  

proc Lookup1D {xVals fxVals x} {
	set idx [BracketLookup $xVals $x]
	set idx1 [lindex $idx 0]; set idx2 [lindex $idx 1]
	set x1 [lindex $xVals $idx1]; set x2 [lindex $xVals $idx2]
	set fx1 [lindex $fxVals $idx1]; set fx2 [lindex $fxVals $idx2]
	set output [Interpolate1D $x1 $x2 $fx1 $fx2 $x]
	return $output
}

# ----------------------------------------------------------------------------------------------------------
# 2D lookup

# Function performs a two dimensional graphical lookup using the input x. Output is bilinearly interpolated.
#
#        | yVals 		
#   -----|-------
#  xVals | fxyVals
#        |
#
# Function takes inputs of xVals and yVals, lists that define the lookup points of the inputs x and y to the function fxy.
# fxyVals is a list where each element is another list containing the lookup values of fxy. fxyVals is indexed against xVals, and the lower
# level lists are indexed against yVals. 
# ** For a full explanation of how to properly use the function and initialise the inputs, refer to the separate file 2DLookupDemo.tcl in the
# same directory as this resource **

proc Lookup2D {xVals yVals fxyVals x y} {
	set xidx [BracketLookup $xVals $x]; set yidx [BracketLookup $yVals $y]
	set xidx1 [lindex $xidx 0]; set xidx2 [lindex $xidx 1]
	set yidx1 [lindex $yidx 0]; set yidx2 [lindex $yidx 1]
	set x1 [lindex $xVals $xidx1]; set x2 [lindex $xVals $xidx2]
	set y1 [lindex $yVals $yidx1]; set y2 [lindex $yVals $yidx2]
	set fx1y1 [lindex [lindex $fxyVals $xidx1] $yidx1]
	set fx1y2 [lindex [lindex $fxyVals $xidx1] $yidx2]
	set fx2y1 [lindex [lindex $fxyVals $xidx2] $yidx1]
	set fx2y2 [lindex [lindex $fxyVals $xidx2] $yidx2]
	set output [Interpolate2D $x1 $y1 $x2 $y2 $fx1y1 $fx1y2 $fx2y1 $fx2y2 $x $y]
	return $output
}

# ----------------------------------------------------------------------------------------------------------  
# Fault integrator Count function

# When trigger is true, the counter counts up for 'confirm' seconds, otherwise it counts down to zero in 'clear' seconds. The countIn
# and prevOutput arguments are the fault integrator count and output from the previous iteration. The function returns the current iteration's 
# count and output as a list, which must be accessed using the lindex command and passed back into the fault integrator function using the countIn
# and prevOutput arguments respectively.

proc FaultInt {trigger confirmTime clearTime countIn prevOutput} {
	global sampTime
	if {$trigger == 1} {
		set countOut [expr {$countIn + $sampTime}]
	} else {
		set countOut [expr {($countIn - [SafeDivide $confirmTime $clearTime]*$sampTime)}]
	}
	set countOut [MyMin $countOut $confirmTime]; set countOut [MyMax $countOut 0]
	set clear [expr {$countOut <= 0}]; set expired [expr {$countOut >= $confirmTime}]
	set output [Latch $expired $clear $prevOutput]
	return [list $output $countOut]
}

# ----------------------------------------------------------------------------------------------------------
# Delay On function

# When trigger is true timer starts counting, when count reaches delay the output goes true. If the trigger goes false the counter resets to zero.

 proc DelayOn {trigger delayTime countIn} {
	global sampTime
	if {$trigger == 1} {
		if {$countIn > $delayTime} {
			set countOut $countIn
		} else {
			set countOut [expr {$countIn + $sampTime}]
		}
		set output [expr {$countOut > $delayTime}]
	} else {
		set countOut 0; set output 0
	}
	return [list $output $countOut]
}

# ----------------------------------------------------------------------------------------------------------
# Delay Off function

# When trigger is false timer starts counting, when count reaches delay the output goes false. If the trigger goes true the counter resets to zero.

 proc DelayOff {trigger holdTime reset countIn} {
	global sampTime
	if { ($reset == 1) || ($trigger == 1)} {
		set output $trigger; set countOut 0
	} else {
		if {$countIn >= $holdTime} {
			set countOut $countIn; set output 0
		} else {
			set countOut [expr {$countIn + $sampTime}]
			if {$countOut >= $holdTime} {
				set output 0
			} else {
				set output 1
			}
		}
	}
	return [list $output $countOut]
}

# ----------------------------------------------------------------------------------------------------------
# False-True Edge Detection

# Function takes inputs of the current and previous values of a boolean signal, if it detects that the signal has transitioned
# from false (0) to true (1) since the previous sample, output is true, otherwise output is False.
 
proc FalseTrueDetect {input prevInput} {
	if {($input == 1) && ($prevInput == 0)} {
		set output 1
	} else {
		set output 0
	}
	return $output
}

# ----------------------------------------------------------------------------------------------------------
# True-False Edge Detection

# Function takes inputs of the current and previous values of a boolean signal, if it detects that the signal has transitioned
# from True (1) to False (0) since the previous sample, output is True, otherwise output is False.

proc TrueFalseDetect {input prevInput} {
	if {($input == 0) && ($prevInput == 1)} {
		set output 1
	} else {
		set output 0
	}
	return $output
}

# ----------------------------------------------------------------------------------------------------------
# Divide by Zero Protection

# Protects against Divide by Zero

proc DividebyZeroProtect {input prevOutput} {
	if {$input == 0} {
		return $prevOutput
	} else {
		return $input
	}
}

# ----------------------------------------------------------------------------------------------------------

namespace export [lsort [info procs [namespace current]::*]]

}