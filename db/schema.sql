CREATE TABLE proteins (
    id BIGSERIAL PRIMARY KEY,
    accession CHAR(10) UNIQUE NOT NULL,
    header TEXT NOT NULL,
    aa_sequence TEXT NOT NULL
);

CREATE INDEX protein_accession_idx ON proteins (accession);


CREATE TABLE peptides (
    id BIGSERIAL NOT NULL,
    aa_sequence CHAR(60) NOT NULL,
    length INTEGER NOT NULL,
    number_of_missed_cleavages SMALLINT NOT NULL,
    weight BIGINT NOT NULL,
    digest_enzym CHAR(5) NOT NULL,
    UNIQUE (aa_sequence, weight),
    PRIMARY KEY (id, weight)
) PARTITION BY RANGE (weight);

-- create partitions
-- CREATE TABLE peptides_%d PARTITION OF peptides FOR VALUES FROM (%d) TO (%d);
CREATE TABLE peptides_0 PARTITION OF peptides FOR VALUES FROM (1) TO (2014015050);
CREATE TABLE peptides_99 PARTITION OF peptides FOR VALUES FROM (2014015051) TO (12000000000);


CREATE INDEX peptide_aa_sequence_idx ON peptides (aa_sequence);
CREATE INDEX peptide_length_idx ON peptides (length);


CREATE TABLE peptides_proteins (
    peptide_id BIGINT NOT NULL,
    protein_id BIGINT NOT NULL,
    PRIMARY KEY (peptide_id, protein_id)
);


CREATE TABLE IF NOT EXISTS amino_acid_modifications (
    id BIGSERIAL PRIMARY KEY,
    accession TEXT NOT NULL,
    name TEXT,
    position SMALLINT NOT NULL,
    is_fix BOOLEAN NOT NULL,
    amino_acid_one_letter_code CHAR(1) NOT null,
    mono_mass BIGINT NOT NULL,
    UNIQUE (accession)
);


CREATE TABLE IF NOT EXISTS base_decoys (
    id BIGSERIAL PRIMARY KEY,
    header TEXT NOT NULL,
    aa_sequence VARCHAR(60) NOT NULL,
    length INTEGER NOT NULL,
    weight BIGINT NOT NULL,
    UNIQUE (aa_sequence, weight)
);


CREATE TABLE IF NOT EXISTS modified_decoys (
    id BIGSERIAL PRIMARY KEY,
    base_decoy_id BIGINT REFERENCES base_decoys(id) NOT NULL,
    c_terminus_modification_id BIGINT REFERENCES amino_acid_modifications(id) NOT NULL,
    n_terminus_modification_id BIGINT REFERENCES amino_acid_modifications(id) NOT NULL,
    modification_ids BIGINT ARRAY NOT NULL,
    weight BIGINT NOT null,
    header_addition TEXT NOT NULL,
    UNIQUE (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight)
);