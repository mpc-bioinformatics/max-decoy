CREATE TABLE proteins (
    id SERIAL PRIMARY KEY,
    accession varchar(10) UNIQUE NOT NULL,
    header text NOT NULL,
    aa_sequence text NOT NULL
);

CREATE INDEX protein_id_idx ON proteins USING gin(id);
CREATE INDEX protein_accession_gin_idx ON proteins USING gin(accession);


CREATE TABLE peptides (
    id SERIAL PRIMARY KEY,
    aa_sequence text UNIQUE NOT NULL,
    length integer NOT NULL,
    number_of_missed_cleavages integer,
    weight FLOAT8 NOT NULL,
    digest_enzym varchar(255) NOT NULL
);

CREATE INDEX peptide_id_idx ON peptides USING gin(id);
CREATE INDEX peptide_weight_gin_idx ON peptides USING gin(weight);
CREATE INDEX peptide_length_gin_idx ON peptides USING gin(length);