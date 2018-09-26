use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;
use proteomic::models::peptide_protein_association::PeptideProteinAssociation;

#[test]
pub fn test_equality() {
    let pep1: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0, 0);
    let pep2: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0, 0);
    let pro1: Protein = Protein::new(
        String::from(">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3"), 
        String::from("MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH")
    );
    let pro2: Protein = Protein::new(
        String::from(">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3"), 
        String::from("MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH")
    );
    let pep_pro_association1 = PeptideProteinAssociation(&pep1, & pro1);
    let pep_pro_association2 = PeptideProteinAssociation(&pep2, & pro2);
    assert!(pep_pro_association1 == pep_pro_association2);
}

#[test]
pub fn test_unequlity() {
    let pep1: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0, 0);
    let pep2: Peptide = Peptide::new(String::from("MVNLLQIVR"), String::from("Trypsin"), 0, 0);
    let pro1: Protein = Protein::new(
        String::from(">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3"), 
        String::from("MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH")
    );
    let pro2: Protein = Protein::new(
        String::from(">sp|O75438|NDUB1_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 1 OS=Homo sapiens OX=9606 GN=NDUFB1 PE=1 SV=1"), 
        String::from("MVNLLQIVRDHWVHVLVPMGFVIGCYLDRKSDERLTAFRNKSMLFKRELQPSEEVTWK")
    );
    let pep_pro_association1 = PeptideProteinAssociation(&pep1, & pro1);
    let pep_pro_association2 = PeptideProteinAssociation(&pep2, & pro2);
    assert!(pep_pro_association1 != pep_pro_association2);
}