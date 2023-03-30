#!/bin/bash

# General settings, customize
PROM="/home/yourusername/bin/promoter-2.0"
SYSTEM="$(uname -s)"

# Select ECHO executable (avoid shell built-in)
if [ "$SYSTEM" == "Linux" ]; then
    ECHO="/bin/echo -e"
else
    ECHO="/bin/echo"
fi

# Parse command line
infile=()
f_opt=""
w_opt=""
V_opt=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        -f)
            f_opt="yes"
            ;;
        -w)
            w_opt="yes"
            ;;
        -V)
            V_opt="yes"
            ;;
        *)
            infile+=("$1")
            ;;
    esac
    shift
done

# Display version info
if [ -n "$V_opt" ]; then
    $ECHO "promoter 2.0a, May 2001, revised in Aug 2004"
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

# Prepare for temporary files
PROMTMP="$PROM/tmp/$$"
mkdir -p "$PROMTMP/data_in" "$PROMTMP/data_out"

# Read input
for f in "${infile[@]}"; do
    {
        echo "$f" | grep -E '^>'
        grep -v -E '^>' "$f"
    } > "$PROMTMP/data_in/$(basename "$f")"
done

cd "$PROM"

# Main
if [ -n "$w_opt" ]; then
    cat "$PROM/html/head.html"
fi

for f in "$PROMTMP/data_in/"*; do
    cat "$f" | "$PROM/bin/promoter_$SYSTEM" >"$PROMTMP/data_out/$(basename "$f").pred"

    if [ -n "$f_opt" ]; then
        $ECHO "INPUT SEQUENCE:"
        cat "$f"

        if [ -n "$w_opt" ]; then
            $ECHO "PREDICTED TRANSCRIPTION START SITES:"
        else
            $ECHO "PREDICTED TRANSCRIPTION START SITES:"
        fi
    fi

    cat "$PROMTMP/data_out/$(basename "$f").pred" | \
    while read -r line; do
        echo "$line"
    done

done

if [ -n "$w_opt" ]; then
    cat "$PROM/html/foot.html"
fi

# Clean up
cd /
rm -r "$PROMTMP"

# End of script
