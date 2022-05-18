create table if not exists genomic_islands
-- Genomic islands
(
    id   integer not null primary key autoincrement,
    name text    not null unique, -- GI name, suffixed with (cGI#) in the event that two different GIs are given the same name in their publications
    role text                     -- Description of the GI and its function
);

create table if not exists alternate_names
-- Alternate island names
(
    gi   integer not null,
    name text    not null,
    foreign key (gi) references genomic_islands (id),
    unique (gi, name)
);

create table if not exists gi_type
-- Genomic Island Types
(
    gi integer not null,
    type text  not null, -- GI type. One of "prophage","ICE","transposon","integron","integrated_plasmid"
    foreign key (gi) references genomic_islands (id),
    unique (gi, type),
    check ( type IN ("prophage","ICE","transposon","integron","integrated_plasmid") )
);

create table if not exists strains
-- Microbe strains
(
    id    integer not null primary key autoincrement,
    gbuid text    unique, -- Genbank unique id for strain
    name  text    not null
);
-- https://stackoverflow.com/a/49846452/15446750
create unique index if not exists `strain_name` on `strains` (
    coalesce(gbuid, name),
    name
);

create table if not exists sequences
-- Genomic island DNA sequences. Multiple variant sequences can be attributed to a single island
(
    id     integer not null primary key autoincrement,
    gi     integer not null,
    gbuid  text unique,             -- Genbank unique id for GI sequence
    start  integer,                 -- File sequence start relative to Genbank hosted sequence
    end    integer,                 -- File sequence end relative to Genbank hosted sequence
    gc     real,                    -- % GC composition
    length integer,                 -- Sequence length
    path   text    not null unique, -- path to fasta relative to this db
    foreign key (gi) references genomic_islands (id)
);

create table if not exists sources
-- Source genomic island was isolated from. A GI can have multiple sources.
(
    id     integer not null primary key autoincrement,
    gi     integer not null,
    strain integer,
    start  integer, -- GI start relative to strain reference genome
    end    integer, -- GI end relative to strain reference genome
    size   integer, -- GI size
    seq    integer unique,
    foreign key (gi) references genomic_islands (id),
    foreign key (strain) references strains (id),
    foreign key (seq) references sequences (id)
);
create index if not exists sources_seq on sources (seq);

create table if not exists publications
-- Publications for sources
(
    id          integer not null primary key autoincrement,
    publication text    not null unique,
    doi         text unique
);

create table if not exists source_pub_assoc
-- Source - Publication associations
(
    source      integer not null,
    publication integer not null,
    foreign key (source) references sources (id),
    foreign key (publication) references publications (id)
);
