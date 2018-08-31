CREATE TABLE proteins (
    id SERIAL PRIMARY KEY,
    accession varchar(10) UNIQUE NOT NULL,
    header text NOT NULL,
    aa_sequence text NOT NULL
);

CREATE INDEX protein_id_idx ON proteins (id);
CREATE INDEX protein_accession_idx ON proteins (accession);


CREATE TABLE peptides (
    id SERIAL PRIMARY KEY,
    aa_sequence varchar(60) UNIQUE NOT NULL,
    length integer NOT NULL,
    number_of_missed_cleavages integer,
    weight integer NOT NULL,
    digest_enzym varchar(60) NOT NULL
);

CREATE INDEX peptide_id_idx ON peptides (id);
CREATE INDEX peptide_weight_idx ON peptides (weight);
CREATE INDEX peptide_length_idx ON peptides (length);