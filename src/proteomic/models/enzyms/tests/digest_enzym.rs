use proteomic::models::persistable::{Persistable};
use proteomic::models::peptide_protein_association::PeptideProteinAssociation;
use proteomic::models::peptides::peptide::Peptide;
use proteomic::models::protein::Protein;
use proteomic::models::enzyms;
use proteomic::models::enzyms::digest_enzym::DigestEnzym;
use proteomic::models::enzyms::digest_summary::DigestSummary;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::amino_acids::amino_acid::AminoAcid;


/// Test protein P77377
const P77377_HEADER: &str = ">sp|P77377|WZXC_ECOLI Lipopolysaccharide biosynthesis protein WzxC OS=Escherichia coli (strain K12) OX=83333 GN=wzxC PE=1 SV=1";
const P77377_SEQUENCE: &str = "MSLREKTISGAKWSAIATVIIIGLGLVQMTVLARIIDNHQFGLLTVSLVIIALADTLSDFGIANSIIQRKEISHLELTTLYWLNVGLGIVVCVAVFLLSDLIGDVLNNPDLAPLIKTLSLAFVVIPHGQQFRALMQKELEFNKIGMIETSAVLAGFTCTVVSAHFWPLAMTAILGYLVNSAVRTLLFGYFGRKIYRPGLHFSLASVAPNLRFGAWLTADSIINYLNTNLSTLVLARILGAGVAGGYNLAYNVAVVPPMKLNPIITRVLFPAFAKIQDDTEKLRVNFYKLLSVVGIINFPALLGLMVVSNNFVPLVFGEKWNSIIPVLQLLCVVGLLRSVGNPIGSLLMAKARVDISFKFNVFKTFLFIPAIVIGGQMAGAIGVTLGFLLVQIINTILSYFVMIKPVLGSSYRQYILSLWLPFYLSLPTLVVSYALGIVLKGQLALGMLLAVQIATGVLAFVVMIVLSRHPLVVEVKRQFCRSEKMKMLLRAG";

/// Resulting peptides for P77377 when digested with Trypsin and digest parameters:
/// * missed cleavages = 2
/// * minimum peptide length = 6
/// * maximum peptide length = 50
const RESULTING_PEPTIDES_FOR_TRYPSIN: [&str; 71] = [
    "IGMIETSAVLAGFTCTVVSAHFWPLAMTAILGYLVNSAVRTLLFGYFGRK",
    "IGMIETSAVLAGFTCTVVSAHFWPLAMTAILGYLVNSAVRTLLFGYFGR",
    "LLSVVGIINFPALLGLMVVSNNFVPLVFGEKWNSIIPVLQLLCVVGLLR",
    "TFLFIPAIVIGGQMAGAIGVTLGFLLVQIINTILSYFVMIKPVLGSSYR",
    "FGAWLTADSIINYLNTNLSTLVLARILGAGVAGGYNLAYNVAVVPPMK",
    "KEISHLELTTLYWLNVGLGIVVCVAVFLLSDLIGDVLNNPDLAPLIK",
    "EISHLELTTLYWLNVGLGIVVCVAVFLLSDLIGDVLNNPDLAPLIK",
    "ELEFNKIGMIETSAVLAGFTCTVVSAHFWPLAMTAILGYLVNSAVR",
    "KIYRPGLHFSLASVAPNLRFGAWLTADSIINYLNTNLSTLVLAR",
    "IYRPGLHFSLASVAPNLRFGAWLTADSIINYLNTNLSTLVLAR",
    "IGMIETSAVLAGFTCTVVSAHFWPLAMTAILGYLVNSAVR",
    "LRVNFYKLLSVVGIINFPALLGLMVVSNNFVPLVFGEK",
    "ILGAGVAGGYNLAYNVAVVPPMKLNPIITRVLFPAFAK",
    "GQLALGMLLAVQIATGVLAFVVMIVLSRHPLVVEVKR",
    "IIDNHQFGLLTVSLVIIALADTLSDFGIANSIIQRK",
    "VNFYKLLSVVGIINFPALLGLMVVSNNFVPLVFGEK",
    "GQLALGMLLAVQIATGVLAFVVMIVLSRHPLVVEVK",
    "IIDNHQFGLLTVSLVIIALADTLSDFGIANSIIQR",
    "WNSIIPVLQLLCVVGLLRSVGNPIGSLLMAKAR",
    "WNSIIPVLQLLCVVGLLRSVGNPIGSLLMAK",
    "LLSVVGIINFPALLGLMVVSNNFVPLVFGEK",
    "ILGAGVAGGYNLAYNVAVVPPMKLNPIITR",
    "EKTISGAKWSAIATVIIIGLGLVQMTVLAR",
    "QYILSLWLPFYLSLPTLVVSYALGIVLK",
    "GQLALGMLLAVQIATGVLAFVVMIVLSR",
    "TLLFGYFGRKIYRPGLHFSLASVAPNLR",
    "TISGAKWSAIATVIIIGLGLVQMTVLAR",
    "TLSLAFVVIPHGQQFRALMQKELEFNK",
    "FGAWLTADSIINYLNTNLSTLVLAR",
    "ILGAGVAGGYNLAYNVAVVPPMK",
    "LNPIITRVLFPAFAKIQDDTEK",
    "WSAIATVIIIGLGLVQMTVLAR",
    "SVGNPIGSLLMAKARVDISFK",
    "TLSLAFVVIPHGQQFRALMQK",
    "KIYRPGLHFSLASVAPNLR",
    "WNSIIPVLQLLCVVGLLR",
    "IYRPGLHFSLASVAPNLR",
    "VLFPAFAKIQDDTEKLR",
    "TLSLAFVVIPHGQQFR",
    "VLFPAFAKIQDDTEK",
    "LNPIITRVLFPAFAK",
    "SVGNPIGSLLMAKAR",
    "IQDDTEKLRVNFYK",
    "HPLVVEVKRQFCR",
    "ARVDISFKFNVFK",
    "SVGNPIGSLLMAK",
    "MSLREKTISGAK",
    "VDISFKFNVFK",
    "ALMQKELEFNK",
    "TLLFGYFGRK",
    "QFCRSEKMK",
    "TLLFGYFGR",
    "HPLVVEVKR",
    "IQDDTEKLR",
    "SEKMKMLLR",
    "HPLVVEVK",
    "RQFCRSEK",
    "ARVDISFK",
    "MKMLLRAG",
    "EKTISGAK",
    "VLFPAFAK",
    "QFCRSEK",
    "LNPIITR",
    "IQDDTEK",
    "LRVNFYK",
    "VDISFK",
    "TISGAK",
    "MSLREK",
    "MKMLLR",
    "ELEFNK",
    "MLLRAG"
];

#[test]
/// Tests digestion of Protein with Trypin. Test protein is P77377 from E. Coli.
/// Digestions parameters:
/// * missed cleavages = 2
/// * minimum peptide length = 6
/// * maximum peptide length = 50
///
/// Test procedure:
///
/// 1. Delete all proteins, peptides and peptide/protein-associations.
/// 2. Digest protein
/// 3. Check if protein was created
/// 4. Check if number of created proteins, reported in digestion summary, is equals actual peptide count of database and length of RESULTING_PEPTIDES_FOR_TRYPSIN (check against RESULTING_PEPTIDES_FOR_TRYPSIN length is only valid for this protein with this digest parameters (manually checked), because peptides can be equal other peptides after geralization).
/// 5. Check if every peptide sequence (gerneralized) is found in database
fn test_digestion_with_trypsin() {
    let conn: postgres::Connection = DatabaseConnection::get_database_connection();
    match Protein::delete_all(&conn) {
        Ok(_) => (),
        Err(err) => panic!("proteomic::models::enzyms::tests::digest_enzym.test_digestion_with_trypsin(): Could not delete all proteins, reason: {}", err)
    };
    match Peptide::delete_all(&conn) {
        Ok(_) => (),
        Err(err) => panic!("proteomic::models::enzyms::tests::digest_enzym.test_digestion_with_trypsin(): Could not delete all peptides, reason: {}", err)
    };
    match PeptideProteinAssociation::delete_all(&conn) {
        Ok(_) => (),
        Err(err) => panic!("proteomic::models::enzyms::tests::digest_enzym.test_digestion_with_trypsin(): Could not delete all peptide/protein-associations, reason: {}", err)
    };
    let mut protein: Protein = Protein::new(P77377_HEADER, P77377_SEQUENCE);
    let mut enzym = enzyms::get("trypsin", &conn, 2, 6, 50);
    let summary: DigestSummary = enzym.digest(&mut protein, 100);
    assert!(summary.has_created_protein());
    let peptide_count = match Peptide::count(&conn) {
        Ok(count) => count,
        Err(query_err) => panic!("proteomic::models::enzyms::tests::digest_enzym.test_digestion_with_trypsin(): Could not get peptide count from database, reason: {}", query_err)
    };
    assert_eq!(peptide_count as usize, summary.get_number_of_created_peptides());
    assert_eq!(summary.get_number_of_created_peptides(), RESULTING_PEPTIDES_FOR_TRYPSIN.len());
    for peptide_seq in RESULTING_PEPTIDES_FOR_TRYPSIN.iter() {
        let generalized_peptide_seq: String = AminoAcid::gerneralize_sequence(*peptide_seq);
        match Peptide::exists_where(&conn, "aa_sequence = $1", &[&generalized_peptide_seq]) {
            Ok(_) => (),
            Err(query_err) => panic!("proteomic::models::enzyms::tests::digest_enzym.test_digestion_with_trypsin(): Could not find peptide '{}' ('{}'), reason: {}", peptide_seq, generalized_peptide_seq, query_err)
        };
    }
}