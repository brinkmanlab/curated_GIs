create table genomic_islands
-- Genomic islands
(
    id integer not null primary key autoincrement,
    name text not null,
    type text,
    role text
);

create table strains
-- Microbe strains
(
    id integer not null primary key autoincrement,
    gbuid text,  -- Genbank unique id for strain
    name text not null
);

create table sequences
-- Genomic island DNA sequences. Multiple variant sequences can be attributed to a single island
(
    id integer not null primary key autoincrement,
    gi integer not null foreign key references genomic_islands (id),
    gbuid text,  -- Genbank unique id for GI sequence
    gc real,  -- % GC composition
    path text not null  -- path to fasta relative to this db
);

create table sources
--
(
    id integer not null primary key autoincrement,
    gi integer not null foreign key references genomic_islands (id),
    strain integer foreign key references strains (id),
    start integer,  -- GI start relative to strain reference genome
    end integer,  -- GI end relative to strain reference genome
    size integer,  -- GI size
    seq integer foreign key references sequences (id),
    pmid integer,  -- PubMed ID
    publication text  -- Citation
);