BEGIN {
	FS=" "; OFS="\t"
}
# qseqid qpid qstart qstop spid cog
pqseqid == $1 && $6 != pcog && $3-pend <= 14000 {
	if (pcog < $6) {
		pair = pcog "_" $6
	} else {
		pair = $6 "_" pcog
	}
  print candidate, pair
}
pqseqid != $1 || $3-pend > 14000 {
	candidate = $1 ":" $3
}
{
  pqseqid=$1; pcog=$6; pstart=$3; pend=$4;
}