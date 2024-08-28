# regex for variable names
NAME="[A-z][A-z0-9_]*"
# regex for expressions
EXPR='\[expr *([A-z0-9_\(\)-+=/ .\*]*)\]'
RE_LOGI='[=!><]'
#EXPR='\[expr ([A-z0-9-+/ \(\)\.\*]*)\]'

# use perl
SEDE="perl -pe"

OPS_PREFIX="model"
declare -a OPS_CMDS=(
  node
  fix
  rigidLink
  element
  section
  patch
  layer
  analyze
  algorithm
  system
)

cat - <<EOF
from math import cos,sin,sqrt,pi
import opensees.openseespy as ops
EOF

cat $1 | {
  $SEDE "s|\t|       |g"
} | {
  $SEDE "s/source LibUnits\.tcl/from opensees.units.english import * # CONVERT-COMPLETE/g"
} | {
  #
  # Procedures
  #
  $SEDE "s/^proc ([A-z0-9_]*) \{(.*)\} \{/def \1(\2): # CONVERT-COMPLETE/g"
} | {
  $SEDE "s|lappend ($NAME) (.*)|\1.append(\2)|g"
} | {
  $SEDE 's|set *([A-z][A-z0-9_]*) ([A-z][A-z.0-9_-]*)$|\1 = '\2'|g if !/#/ && !/"/'
} | {
  $SEDE "s|set *($NAME) (.*)|\1 = \2|g"' if !/#/ && !/"/'
} | {
  $SEDE "s|\[incr ($NAME)\]|\1 += 1|g"
} | {
  # when `incr` found inside braces, keep braces
  # so for loops can be identified later
  $SEDE "s|{incr ($NAME) ([A-z0-9_]*)}|{\1 += \2}|g"
} | {
  $SEDE "s|{incr ($NAME) ([A-z0-9_]*)}|{\1 += \2}|g"
} | {
#
## wrap words in quotes if line does not contain a " or start with #
#
  # end of line
  $SEDE "s:([A-z'])"' ([A-z][0-9A-z_/\.]{1,})$'":\1 '\2':g if "' !/"/ && !/^ *#/ && !/CONVERT-COMPLETE/ && !/^ *set /'
} | {
  # inside line
  #$SEDE "s:([!\[][A-z']*)"' ([A-z][0-9A-z_/\.]{1,}) '":\1 '\2', :g if /^[ A-z]/ &&"' !/"/ && !/^ *#/ && !/CONVERT-COMPLETE/'
  $SEDE "s:([A-z'])"' ([A-z][0-9A-z_/\.]{1,}) '":\1 '\2', :g if !/\[\]/ &&"' !/"/ && !/^ *#/ && !/CONVERT-COMPLETE/ && !/^ *set /'
} | {

## remove dollar signs
  sed 's:\$::g'
} | {
  # change || to `or`
  sed 's: || : or :g'
} | {
  $SEDE 's:pow\(([A-z0-9.]*),([A-z0-9.]*)\):\1**\2:g'
} | {
  $SEDE "s:set *($NAME) *$EXPR:\1 = \2:g"
} | {
  $SEDE 's:( *)set *([A-z0-9_]*) *\[expr ([^]]+)\]:$1$2 = $3:g'
} | {
  $SEDE "s:set *($NAME) *([-0-9.]*):\1 = \2:g"
} | {
  $SEDE 's: \[expr ([^]]*)\];:, $1:g'
} | {
  $SEDE 's: \[expr ([^]]*)\] :, ($1) :g'
} | {
  $SEDE 's: \[expr ([-+/* A-z0-9_()]+)\]:, ($1):g'
} | {
#
# Control flow
#
  # if
  $SEDE 's|^([ \t]*)if {([A-z0-9-+"!=/ \.\*]*)} {$|\1if \2:|g'
} | {
  # else
  $SEDE "s|} else {|else:|g"
} | {
  # for {i = start} {test} {incr}
  $SEDE "s|for \{($NAME) *= *0\} \{$NAME *$RE_LOGI *($NAME)\} \{$NAME *+= *1\} \{|for \1 in range(\2):|g"
} | {
  $SEDE "s|for {i = 0} {i < 8} {i += 1}|for i in range(8):|g"
} | {
  # remove closing brace
  $SEDE "s|^([ \t]*)};* *$||g"
} | {
#
# Misc
#
  # && Operator
  $SEDE 's|&&| and |g'
} | {
  # puts command with form `puts "string"`
  $SEDE 's|puts "(.*)"|print("\1")|g'

  # replace use of *in unit
  $SEDE 's|\*in([/\*-+ ])|*inch\1|g'
} | {
#
# Flag options
#
  # first, treat flag followed by flag as boolean
  $SEDE s/" -([A-z]*) -"/' \1=True, -/g if !/"/ && !/#/'
} | {
  # treat flags at end of line as boolean
  $SEDE s/" -([A-z]*)$"/' \1=True/g if !/"/ && !/#/'
} | {
  $SEDE s/" -([A-z]*) ([0-9] [0-9 ]{1,})"/' \1=[\2], /g if !/"/ && !/#/ && !/\(/'
} | {
  $SEDE s/" -([A-z]*) '/ \1='\2/g"' if !/"/ && !/#/ && !/\(/'
} | {
  $SEDE s/" -([A-z]*) ([A-z0-9.]*)"/' \1=\2, /g if !/"/ && !/#/ && !/\(/'
} | {
  $SEDE s/",  section ([0-9]* [ 0-9]*)/, section=[\1], /g"
} | {
  # remove comments at end of line
  $SEDE s/"  *#.*$"//g
} | {
#
# 
# Simple model commands
$SEDE "s/(section) '([A-z0-9]*)', ([A-z 0-9+-.'=,]*)\{"'/${1}.${2}(${3}, [/g if !/"/ && !/#/'
} | {
  $SEDE 's/\[;$/[/g'
} | {
  # commands that probably take the rest of their line
  $SEDE "s/(node|mass|geomTransf|element|uniaxialMaterial) ([A-z \(\)0-9-.'=,*\/]*)"'/model.\1(${2})/g if !/"/ && !/^ *#/'
} | {
  # commands that probably dont have parenthesized expressions
  $SEDE "s/(fix|fixZ|fixY|fixX) ([A-z 0-9-.'=,]*);*"'/model.\1(${2})/g if !/"/ && !/^ *#/'
} | {
  # Simple analysis commands
  $SEDE "s/(system|numberer|recorder|rayleigh|analyze|loadConst|algorithm|test|timeSeries|integrator|analysis|constraints|remove) ([A-z 0-9-.'=,]*)"'/model.\1(${2})/g if !/"/ && !/#/'
} | {
  $SEDE "s/(patch|layer|fiber) *'([A-z]*)', *([A-z 0-9-.'=,]*);{,1}"'/\1.${2}(${3}),/g if !/"/'
} | {
# spaces to commas
  $SEDE "s/([A-z0-9'.]) {1,}([A-z0-9'.])/\1, \2/g"' if /^[A-z]/ && !/"/ && !/^ *#/ && !/^if / && !/CONVERT-COMPLETE/'
} | {
  $SEDE "s/(?=^(def) )([A-z0-9'.]) {1,}([A-z0-9'.])/\1, \2/g"
} | {
  $SEDE "s/=,/=/g"
} 

