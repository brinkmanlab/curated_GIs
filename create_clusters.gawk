BEGIN {
	FS=" "; OFS="\t"
}
# qpid qseqid qstart qstop spid cog
pqseqid == $2 && $6 != pcog && $3-pend <= 14000 {
	if (pcog < $6) {
		pair = pcog "_" $6
	} else {
		pair = $6 "_" pcog
	}
  print $2, pair
}
{
  pqseqid=$2; pcog=$6; pstart=$3; pend=$4;
}