# qpid qseqid qstart qstop spid? cog?

function print_gi(seqid, start, end, tokens) {
	label=seqid ":" start "-" end;
	print label, tokens[0];
	for (i=1; i < length(tokens); i++) {
	  if (tokens[i-1] < tokens[i]) {
	  	print label, tokens[i-1] "_" tokens[i];
	  } else {
	  	print label, tokens[i] "_" tokens[i-1];
	  }
	  print label, tokens[i];
	}
}

function trim_left(tokens, until) {
	len = length(tokens);
	for (i=0; i + until < len; i++) {
		tokens[i] = tokens[i + until];
	}
	for (i=len-1; len - i <= until; i--) {
		delete tokens[i];
	}
}

BEGIN {
	FS=" "; OFS="\t";
	PROCINFO["sorted_in"] = "@ind_num_asc";
	delete tokens[0];
	insertion_coeff = (1/cutoff)-1;  # There are a maximum of (#pcogs * insertion_coeff) insertions allowed within a GI
}

# GI starts at first blast hit
pcog && !GI_seqid {
	GI_seqid = pqseqid;
	GI_start = pstart;
	insertions = 0;
}

# Skip duplicate record, retaining the assigned cog if any
pqseqid == GI_seqid &&
pstart == $3 {
	if ($6) pcog = $6;
	next;
}

GI_seqid && pqseqid == GI_seqid {
  if (pcog) {
  	i = length(tokens);
  	tokens[i] = pcog;
  	lastpcog[length(lastpcog)] = i;
  	lastend = pend;
  } else {
    insertions++;
    tokens[i] = pqpid;
	}
	while (length(lastpcog) * insertions_coeff < insertions) {
		# reduce insertions
		# output candidate gi bounded by all tokens
		for (i = 0; i <= lastpcog[length(lastpcog)-1]; ++i) candidate[i] = tokens[i];
		print_gi(GI_seqid, GI_start, lastend, candidate);
		trim_left(tokens, lastpcog[0]);
		trim_left(lastpcog, 1);
	}
}

# End of GI encountered, output tokens
($2 != GI_seqid || $3-pend > 14000) && GI_seqid {
	print_gi(GI_seqid, GI_start, pend, tokens);
	delete tokens;
	GI_seqid = "";
}

{
  pqpid=$1;
  pqseqid=$2;
  pstart=$3;
  pend=$4;
	pcog=$6;
}

END {
  if (GI_seqid) {
		tokens[length(tokens)] = pcog ? pcog : pqpid;
		print_gi(GI_seqid, GI_start, pend, tokens);
	}
}