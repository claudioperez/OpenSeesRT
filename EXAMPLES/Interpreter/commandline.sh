
# run file, then drop into a REPL
python -m opensees -i truss.tcl

# run file 'args.tcl'. Output should be 'a b c'
python -m opensees args.tcl a b c

# This should be a file not found error
cat args.tcl | python -m opensees a b c

# run command 'set argv 3'. Output should be '3' on stdout
python -m opensees -c 'puts "[set argv 3]"'

# run command 'set argv 3', then read file 'args.tcl', then drop into a REPL.
# The command should override argv and the output should be '3' on stdout
python -m opensees -c 'set argv 3' args.tcl a b c


