use proteomic::models::protein::Protein;

#[test]
pub fn test_equality() {
    let pro1: Protein = Protein::new(
        ">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3",
        "MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH"
    );
    let pro2: Protein = Protein::new(
        ">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3",
        "MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH"
    );
    assert!(pro1 == pro2);
}

#[test]
pub fn test_unequlity() {
    let pro1: Protein = Protein::new(
        ">sp|O95139|NDUB6_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 6 OS=Homo sapiens OX=9606 GN=NDUFB6 PE=1 SV=3",
        "MTGYTPDEKLRLQQLRELRRRWLKDQELSPREPVLPPQKMGPMEKFWNKFLENKSPWRKMVHGVYKKSIFVFTHVLVPVWIIHYYMKYHVSEKPYGIVEKKSRIFPGDTILETGEVIPPMKEFPDQHH"
    );
    let pro2: Protein = Protein::new(
        ">sp|O75438|NDUB1_HUMAN NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 1 OS=Homo sapiens OX=9606 GN=NDUFB1 PE=1 SV=1",
        "MVNLLQIVRDHWVHVLVPMGFVIGCYLDRKSDERLTAFRNKSMLFKRELQPSEEVTWK"
    );
    assert!(pro1 != pro2);
}