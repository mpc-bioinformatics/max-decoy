CREATE TABLE proteins (
    id SERIAL PRIMARY KEY,
    accession char(10) UNIQUE NOT NULL,
    header text NOT NULL,
    aa_sequence text NOT NULL
);

CREATE INDEX protein_accession_idx ON proteins (accession);


CREATE TABLE peptides (
    id SERIAL PRIMARY KEY,
    aa_sequence char(60) UNIQUE NOT NULL,
    length integer NOT NULL,
    number_of_missed_cleavages integer,
    weight integer NOT NULL,
    digest_enzym char(60) NOT NULL
);

CREATE INDEX peptide_aa_sequence_idx ON peptides (aa_sequence);
CREATE INDEX peptide_weight_idx ON peptides (weight);
CREATE INDEX peptide_length_idx ON peptides (length);

CREATE TABLE peptides_proteins (
    peptide_id INTEGER,
    protein_id INTEGER,
    PRIMARY KEY (peptide_id, protein_id)
);


CREATE OR REPLACE FUNCTION self_wipe()
    RETURNS void
AS $$
    DECLARE
    BEGIN
        DROP TABLE IF EXISTS proteins;
        DROP TABLE IF EXISTS peptides;
        DROP TABLE IF EXISTS peptides_proteins;

        CREATE TABLE proteins (
            id SERIAL PRIMARY KEY,
            accession char(10) UNIQUE NOT NULL,
            header text NOT NULL,
            aa_sequence text NOT NULL
        );

        CREATE INDEX protein_accession_idx ON proteins (accession);

        CREATE TABLE peptides (
            id SERIAL PRIMARY KEY,
            aa_sequence char(60) UNIQUE NOT NULL,
            length integer NOT NULL,
            number_of_missed_cleavages integer,
            weight integer NOT NULL,
            digest_enzym char(60) NOT NULL
        );

        CREATE INDEX peptide_weight_idx ON peptides (weight);
        CREATE INDEX peptide_length_idx ON peptides (length);

        CREATE TABLE peptides_proteins (
            peptide_id INTEGER,
            protein_id INTEGER,
            PRIMARY KEY (peptide_id, protein_id)
        );

    END;
$$ LANGUAGE plpgsql;