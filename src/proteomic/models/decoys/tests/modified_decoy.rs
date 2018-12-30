use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::persistable::Persistable;
use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::decoys::modified_decoy::ModifiedDecoy;
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::amino_acids::modification::Modification;
use proteomic::models::amino_acids::modification::ModificationPosition;

#[test]
fn modified_decoy_insert() {
    let conn = DatabaseConnection::get_database_connection();
    ModifiedDecoy::delete_all(&conn);
    BaseDecoy::delete_all(&conn);
    Modification::delete_all(&conn);
    assert_eq!(ModifiedDecoy::get_count(&conn), 0);
    let aa_sequence: &str = "CGKMM";
    let mut base_decoy: BaseDecoy = BaseDecoy::new(
        ">DECOY|ModifiedDecoy insert test",
        aa_sequence,
        AminoAcid::get_sequence_weight(aa_sequence)
    );
    match base_decoy.create(&conn) {
        Ok(_) => println!("BaseDecoy successfully created"),
        Err(err) => panic!("BaseDecoy not created, original error {}", err)
    }
    let mut modification: Modification = Modification::new(
        "modified_decoy_test:insert", 
        "ModifiedDecoy inert test",
        ModificationPosition::Anywhere,
        true,
        'C',
        25.1
    );
    match modification.create(&conn) {
        Ok(_) => println!("Modification successfully created"),
        Err(err) => panic!(err)
    }
    let mut modified_decoy: ModifiedDecoy = ModifiedDecoy::new(&base_decoy);
    match modified_decoy.set_modification_at(0, Some(modification.clone())) {
        Ok(_) => println!("Modification applied"),
        Err(err) => panic!("Could not set Modification on ModifiedDecoy, original error: {}", err)
    }
    assert_eq!(base_decoy.get_weight() + modification.get_mono_mass(), modified_decoy.get_weight());
    match modified_decoy.create(&conn) {
        Ok(_) => println!("ModifiedDecoy successfully created"),
        Err(err) => panic!("ModifiedDecoy not created, original error {}", err)
    }
    assert!(modified_decoy.get_primary_key() > 0);
    assert_eq!(ModifiedDecoy::get_count(&conn), 1);
}
