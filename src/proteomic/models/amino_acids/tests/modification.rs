use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::persistable::Persistable;
use proteomic::models::amino_acids::modification::Modification;
use proteomic::models::amino_acids::modification::ModificationPosition;

#[test]
fn modification_insert() {
    let conn = DatabaseConnection::get_database_connection();
    Modification::delete_all(&conn);
    let mut modification_count: i64 = match Modification::count(&conn) {
        Ok(count) => count,
        Err(err) => panic!("proteomic::models::amino_acids::tests::modification::modification_insert at first Modification::count(): {}", err)
    };
    assert_eq!(modification_count, 0);
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
        Err(err) => panic!("proteomic::models::amino_acids::tests::modification::modification_insert at Modification.create(): {}", err)
    }
    assert!(modification.get_primary_key() > 0);
    modification_count = match Modification::count(&conn) {
        Ok(count) => count,
        Err(err) => panic!("proteomic::models::amino_acids::tests::modification::modification_insert at first Modification::get_count(): {}", err)
    };
    assert_eq!(modification_count, 1);
}