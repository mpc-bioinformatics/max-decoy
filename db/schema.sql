CREATE TABLE proteins (
    id SERIAL PRIMARY KEY,
    accession varchar(10) UNIQUE NOT NULL,
    header text NOT NULL,
    aa_sequence text NOT NULL
);

CREATE INDEX protein_accession_idx ON proteins (accession);


CREATE TABLE peptides (
    id SERIAL PRIMARY KEY,
    aa_sequence varchar(60) UNIQUE NOT NULL,
    length integer NOT NULL,
    number_of_missed_cleavages integer,
    weight integer NOT NULL,
    digest_enzym varchar(60) NOT NULL
);

CREATE INDEX peptide_weight_idx ON peptides (weight);
CREATE INDEX peptide_length_idx ON peptides (length);

CREATE TABLE peptides_proteins (
    peptide_aa_sequence varchar(60),
    protein_accession varchar(10),
    PRIMARY KEY (peptide_aa_sequence, protein_accession)

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
            accession varchar(10) UNIQUE NOT NULL,
            header text NOT NULL,
            aa_sequence text NOT NULL
        );

        CREATE INDEX protein_accession_idx ON proteins (accession);

        CREATE TABLE peptides (
            id SERIAL PRIMARY KEY,
            aa_sequence varchar(60) UNIQUE NOT NULL,
            length integer NOT NULL,
            number_of_missed_cleavages integer,
            weight integer NOT NULL,
            digest_enzym varchar(60) NOT NULL
        );

        CREATE INDEX peptide_weight_idx ON peptides (weight);
        CREATE INDEX peptide_length_idx ON peptides (length);

        CREATE TABLE peptides_proteins (
            peptide_aa_sequence varchar(60),
            protein_accession varchar(10),
            PRIMARY KEY (peptide_aa_sequence, protein_accession)

        );
    END;
$$ LANGUAGE plpgsql;