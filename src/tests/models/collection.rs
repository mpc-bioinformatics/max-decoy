use proteomic::models::collection::Collection;
use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;

#[test]
pub fn test_len_with_peptides() {
    let mut peptide_collection: Collection<Peptide> = Collection::new();
    let pep1: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0, 0);
    let pep2: Peptide = Peptide::new(String::from("MVNLLQIVR"), String::from("Trypsin"), 0, 0);   
    let pep3: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0, 0);
    peptide_collection.add(pep1);
    assert_eq!(peptide_collection.len(), 1);
    peptide_collection.add(pep2);
    assert_eq!(peptide_collection.len(), 2);
    peptide_collection.add(pep3);
    assert_eq!(peptide_collection.len(), 2);
}

#[test]
pub fn test_len_with_proteins() {
    let mut peptide_collection: Collection<Protein> = Collection::new();
    let pro1: Protein = Protein::new(
        String::from(">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3"), 
        String::from("MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH")
    );
    let pro2: Protein = Protein::new(
        String::from(">sp|O75438|NDUB1_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 1 OS=Homo sapiens OX=9606 GN=NDUFB1 PE=1 SV=1"), 
        String::from("MVNLLQIVRDHWVHVLVPMGFVIGCYLDRKSDERLTAFRNKSMLFKRELQPSEEVTWK")
    );
    let pro3: Protein = Protein::new(
        String::from(">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3"), 
        String::from("MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH")
    );
    peptide_collection.add(pro1);
    assert_eq!(peptide_collection.len(), 1);
    peptide_collection.add(pro2);
    assert_eq!(peptide_collection.len(), 2);
    peptide_collection.add(pro3);
    assert_eq!(peptide_collection.len(), 2);
}