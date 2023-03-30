#!/bin/bash

###############################################################################
# GENERAL SETTINGS, CUSTOMIZE
###############################################################################

# Full path to promoter directory
PROM="/home/user/bin/promoter-2.0"

# Determine platform (do not change this unless you don't have 'uname'!)
SYSTEM=$(uname -s)

# Select AWK executable
case $SYSTEM in
    "IRIX" | "IRIX64" | "AIX" | "OSF1" | "SunOS" )
        AWK="/usr/bin/nawk"
        ;;
    "Linux" )
        AWK="/usr/bin/gawk"
        ;;
    *)
        AWK="mysuperawk"
        ;;
esac

# Select ECHO executable (avoid shell built-in)
if [ $SYSTEM == "Linux" ]; then
    ECHO="/bin/echo -e"
else
    ECHO="/bin/echo"
fi

###############################################################################
# NOTHING SHOULD NEED CHANGING BELOW THIS LINE!
###############################################################################

VER="2.0a"
VERDATE="May 2001, revised in Aug 2004"

# Parse command line
infile=()
f_opt=""
w_opt=""
V_opt=""
for w in "$@"
do
    case $w in
        "-f")
            f_opt="1"
            ;;
        "-w")
            w_opt="1"
            ;;
        "-V")
            V_opt="1"
            ;;
        *)
            infile+=("$w")
            ;;
    esac
done

# Display version info
if [ -n "$V_opt" ]; then
    $ECHO "promoter $VER, $VERDATE"
    exit 0
fi

# Check for binaries
if [ ! -x "$PROM/bin/promoter_$SYSTEM" ]; then
    $ECHO "promoter: cannot find binaries for $SYSTEM"
    exit -1
fi

# Check for input files
for f in "${infile[@]}"; do
    if [ ! -e "$f" ]; then
        $ECHO "promoter: \"$f\" not found"
        exit -2
    fi
done

# Define header flanks
if [ -n "$w_opt" ]; then
    l="<h3>"
    r="</h3>"
else
    l="\n\n"
    r="\n"
fi

if [ -n "$f_opt" ] && [ -n "$w_opt" ]; then
    x="<h3>"
else
    x="$l"
fi

# Prepare for temporary files
PROMTMP="$PROM/tmp/$$"
mkdir -p "$PROMTMP/data_in" "$PROMTMP/data_out"

# Read input
cat "${infile[@]}" | tr '\r' '\n' | grep -v '^$' | \
    $AWK -v DIR
# Read input
cat "${infile[@]}" | tr '\r' '\n' | grep -v '^$' | \
    $AWK -v DIR="$PROMTMP/data" -f "$PROM/bin/fasta2dir"

cd "$PROM"

# Main
if [ -n "$w_opt" ]; then
    cat "$PROM/html/head.html"
fi

for f in "$PROMTMP/data_in/"*; do
    cat "$f" | "$PROM/bin/promoter_$SYSTEM" >"$PROMTMP/data_out/${f##*/}.pred"

    if [ -n "$f_opt" ]; then
        $ECHO "${l}INPUT SEQUENCE:$r"
        cat "$f"

        if [ -n "$w_opt" ]; then
            $ECHO "${l}PREDICTED TRANSCRIPTION START SITES:$r"
        else
            $ECHO "${l}PREDICTED TRANSCRIPTION START SITES:"
        fi
    fi

    cat "$PROMTMP/data_out/${f##*/}.pred" | \
        $AWK -v LEN=$(cat "$PROMTMP/data_out/${f##*/}") -v L="$x" -v R="$r" \
             -v LIN=$(cat "$PROMTMP/data_out/${f##*/}.pred" | wc -l) \
             -f "$PROM/bin/pp"
done

if [ -n "$w_opt" ]; then
    cat "$PROM/html/foot.html"
fi

# Clean up
cd /
rm -r "$PROMTMP"

# End of script
