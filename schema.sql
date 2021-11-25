create table if not exists genomic_islands
-- Genomic islands
(
    id   integer not null primary key autoincrement,
    name text    not null unique, -- GI name, suffixed with (cGI#) in the event that two different GIs are given the same name in their publications
    type text,                    -- TODO CHECK(type IN ("","")  Need to standardise types
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

create table if not exists strains
-- Microbe strains
(
    id    integer not null primary key autoincrement,
    gbuid text unique, -- Genbank unique id for strain
    name  text    not null
);

create table if not exists sequences
-- Genomic island DNA sequences. Multiple variant sequences can be attributed to a single island
(
    id    integer not null primary key autoincrement,
    gi    integer not null,
    gbuid text unique,             -- Genbank unique id for GI sequence
    gc    real,                    -- % GC composition
    path  text    not null unique, -- path to fasta relative to this db
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

create table if not exists pmids
-- PubMed IDs for sources
(
    source integer not null,
    pmid   text    not null, -- PubMed ID
    foreign key (source) references sources (id)
);

create table if not exists publications
-- Publications for sources
(
    source      integer not null,
    publication text    not null,
    foreign key (source) references sources (id)
);