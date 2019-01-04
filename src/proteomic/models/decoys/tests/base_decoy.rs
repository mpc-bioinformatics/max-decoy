use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::persistable::Persistable;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::amino_acids::amino_acid::AminoAcid;

#[test]
fn base_decoy_insert() {
    let conn = DatabaseConnection::get_database_connection();
    BaseDecoy::delete_all(&conn);
    let mut base_decoy_count: i64 = match BaseDecoy::get_count(&conn) {
        Ok(count) => count,
        Err(err) => panic!("{}", err)
    };
    assert_eq!(base_decoy_count, 0);
    let aa_sequence: &str = "CGKMM";
    let mut base_decoy: BaseDecoy = BaseDecoy::new(
        ">DECOY|BaseDecoy insert test",
        aa_sequence,
        AminoAcid::get_sequence_weight(aa_sequence)
    );
    match base_decoy.create(&conn) {
        Ok(_) => println!("BaseDecoy successfully created"),
        Err(err) => panic!("{}", err)
    }
    assert!(base_decoy.get_primary_key() > 0);
    base_decoy_count = match BaseDecoy::get_count(&conn) {
        Ok(count) => count,
        Err(err) => panic!("proteomic::models::decoys::tests::base_decoy::base_decoy_insert at second BaseDecoy::get_count(): {}", err)
    };
    assert_eq!(base_decoy_count, 1);
}

#[test]
fn base_decoy_insert_two_equal() {
    let conn = DatabaseConnection::get_database_connection();
    BaseDecoy::delete_all(&conn);
    let mut base_decoy_count: i64 = match BaseDecoy::get_count(&conn) {
        Ok(count) => count,
        Err(err) => panic!("proteomic::models::decoys::tests::base_decoy::base_decoy_insert at first BaseDecoy::get_count(): {}", err)
    };
    assert_eq!(base_decoy_count, 0);
    let aa_sequence: &str = "CGKMM";
    let mut base_decoy1: BaseDecoy = BaseDecoy::new(
        ">DECOY|BaseDecoy insert test",
        aa_sequence,
        AminoAcid::get_sequence_weight(aa_sequence)
    );
    let mut base_decoy2: BaseDecoy = BaseDecoy::new(
        ">DECOY|BaseDecoy insert test",
        aa_sequence,
        AminoAcid::get_sequence_weight(aa_sequence)
    );
    match base_decoy1.create(&conn) {
        Ok(_) => println!("BaseDecoy successfully created"),
        Err(err) => panic!("proteomic::models::decoys::tests::base_decoy::base_decoy_insert at first BaseDecoy.create(): {}", err)
    }
    assert!(base_decoy1.get_primary_key() > 0);
    match base_decoy2.create(&conn) {
        Ok(_) => println!("BaseDecoy successfully created"),
        Err(err) => panic!("proteomic::models::decoys::tests::base_decoy::base_decoy_insert at first BaseDecoy.create(): {}", err)
    }
    assert!(base_decoy2.get_primary_key() > 0);
    assert_eq!(base_decoy2.get_primary_key(), base_decoy2.get_primary_key());
    base_decoy_count = match BaseDecoy::get_count(&conn) {
        Ok(count) => count,
        Err(err) => panic!("proteomic::models::decoys::tests::base_decoy::base_decoy_insert at second BaseDecoy::get_count(): {}", err)
    };
    assert_eq!(base_decoy_count, 1);
}