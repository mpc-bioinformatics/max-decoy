CREATE TABLE proteins (
    id SERIAL PRIMARY KEY,
    accession CHAR(10) UNIQUE NOT NULL,
    header TEXT NOT NULL,
    aa_sequence TEXT NOT NULL
);

CREATE INDEX protein_accession_idx ON proteins (accession);


CREATE TABLE peptides (
    id SERIAL NOT NULL,
    aa_sequence CHAR(60) NOT NULL,
    length INTEGER NOT NULL,
    number_of_missed_cleavages SMALLINT NOT NULL,
    weight BIGINT NOT NULL,
    digest_enzym CHAR(5) NOT NULL,
    UNIQUE (aa_sequence, weight),
    PRIMARY KEY (id, weight)
) PARTITION BY RANGE (weight);

-- create partitions


CREATE INDEX peptide_aa_sequence_idx ON peptides (aa_sequence);
CREATE INDEX peptide_length_idx ON peptides (length);

CREATE TABLE peptides_proteins (
    peptide_id INTEGER,
    protein_id INTEGER,
    PRIMARY KEY (peptide_id, protein_id)
);
