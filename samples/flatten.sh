# USAGE : sh flatten.sh input.fasta > output_fasta

awk '
/^>/ {
    if (seq) {
        print seq
        seq = ""
    }
    print
}
/^[^>]/ {
    seq = seq $0
}
END {
    if (seq) {
        print seq
    }
}' $1