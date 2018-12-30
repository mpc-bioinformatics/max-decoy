use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::persistable::Persistable;
use proteomic::models::amino_acids::modification::Modification;
use proteomic::models::amino_acids::modification::ModificationPosition;

#[test]
fn modification_insert() {
    let conn = DatabaseConnection::get_database_connection();
    Modification::delete_all(&conn);
    assert_eq!(Modification::get_count(&conn), 0);
    let mut modification: Modification = Modification::new(
        "modification_test:insert", 
        "Modification insert test",
        ModificationPosition::Anywhere,
        true,
        'C',
        25.1
    );
    match modification.create(&conn) {
        Ok(_) => println!("Modification successfully created"),
        Err(err) => panic!(err)
    }
    assert!(modification.get_primary_key() > 0);
    assert_eq!(Modification::get_count(&conn), 1);
}