"""Self-contained calibration set for modality.py.

Run:  python modality.py --selftest

Every SMILES here was retrieved from ChEMBL and its length verified against
`MAX(LENGTH(canonical_smiles))` from an independent query -- see the SMILES
truncation failure mode in SKILL.md. Ground truth is the known chemistry of the
molecule, NOT ChEMBL's `molecule_type`; where the two disagree that is the point
of the case.

The full test is larger (256 cases: 98 IL-17C macrocyclic peptides, 117 IL-17A
compounds, 16 named drugs, 22 oligosaccharides, 3 no-SMILES drugs). It scored
256/256 coarse (small molecule vs not) with 0 false small-molecule calls, and
255/256 fine-grained; the single fine-grained error is ROMIDEPSIN, included
below. This file is the portable subset that travels with the skill.
"""

CASES = [
    # ORITAVANCIN: glycopeptide; ChEMBL `Small molecule`; backbone fraction 0.20.
    ('ORITAVANCIN', 'CN[C@H](CC(C)C)C(=O)N[C@H]1C(=O)N[C@@H](CC(N)=O)C(=O)N[C@H]2C(=O)N[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)O)c4cc(O)cc(O)c4-c4cc3ccc4O)[C@H](O[C@H]3C[C@](C)(N)[C@@H](O)[C@H](C)O3)c3ccc(c(Cl)c3)Oc3cc2cc(c3O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O[C@H]2C[C@](C)(NCc3ccc(-c4ccc(Cl)cc4)cc3)[C@@H](O)[C@H](C)O2)Oc2ccc(cc2Cl)[C@H]1O',
     'Small molecule', 'MOL', 'macrocyclic_peptide'),
    # PLECANATIDE: 
    ('PLECANATIDE', 'CC(C)C[C@H](NC(=O)[C@@H]1CSSC[C@@H]2NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(=O)O)NC(=O)[C@@H](N)CC(N)=O)CSSC[C@H](NC(=O)[C@H](C)NC(=O)[C@H](C(C)C)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](C(C)C)NC2=O)C(=O)N[C@@H]([C@@H](C)O)C(=O)NCC(=O)N1)C(=O)O',
     'Protein', 'BOTH', 'macrocyclic_peptide'),
    # LINACLOTIDE: 
    ('LINACLOTIDE', 'C[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@H](CC(N)=O)NC(=O)[C@@H]2CSSC[C@H](N)C(=O)N[C@H]3CSSC[C@H](NC1=O)C(=O)N[C@@H]([C@@H](C)O)C(=O)NCC(=O)N[C@H](C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)O)CSSC[C@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(=O)O)NC3=O)C(=O)N2',
     'Protein', 'BOTH', 'macrocyclic_peptide'),
    # BORTEZOMIB: HARD NEGATIVE. A dipeptidyl boronate, MW 384 -- 1 alpha linkage. Must NOT be called a peptide.
    ('BORTEZOMIB', 'CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1cnccn1)B(O)O',
     'Small molecule', 'MOL', 'small_molecule'),
    # ROMIDEPSIN: KNOWN GAP. A bicyclic DEPSIpeptide: the backbone alternates amide and ester, so only 2 alpha linkages are found and structure alone says `small_molecule`. It is caught only because ChEMBL types it `Protein`. A depsipeptide typed `Small molecule` would pass as a small molecule.
    ('ROMIDEPSIN', 'C/C=C1\\NC(=O)[C@H]2CSSCC/C=C/[C@H](CC(=O)N[C@H](C(C)C)C(=O)N2)OC(=O)[C@H](C(C)C)NC1=O',
     'Protein', 'MOL', 'macrocyclic_peptide'),
    # MOTIXAFORTIDE: 
    ('MOTIXAFORTIDE', 'N=C(N)NCCC[C@H](NC(=O)[C@@H]1CSSC[C@H](NC(=O)[C@H](Cc2ccc3ccccc3c2)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)c2ccc(F)cc2)C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@H](CCCNC(N)=O)C(=O)N[C@@H](CCCCN)C(=O)N[C@H](CCCCN)C(=O)N2CCC[C@H]2C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CCCNC(N)=O)C(=O)N1)C(N)=O',
     'Protein', 'MOL', 'macrocyclic_peptide'),
    # VOCLOSPORIN: 
    ('VOCLOSPORIN', 'C=C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](C(C)C)C(=O)N1C',
     'Protein', 'MOL', 'macrocyclic_peptide'),
    # DESMOPRESSIN: 
    ('DESMOPRESSIN', 'N=C(N)NCCC[C@@H](NC(=O)[C@@H]1CCCN1C(=O)[C@@H]1CSSCCC(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H](Cc2ccccc2)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CC(N)=O)C(=O)N1)C(=O)NCC(N)=O',
     'Protein', 'MOL', 'macrocyclic_peptide'),
    # DAPTOMYCIN: cyclic lipopeptide; ChEMBL `Small molecule`.
    ('DAPTOMYCIN', 'CCCCCCCCCC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H]1C(=O)NCC(=O)N[C@@H](CCCN)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H](C)C(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)N[C@H](CO)C(=O)N[C@@H]([C@H](C)CC(=O)O)C(=O)N[C@@H](CC(=O)c2ccccc2N)C(=O)O[C@@H]1C',
     'Small molecule', 'MOL', 'macrocyclic_peptide'),
    # TELAVANCIN: 
    ('TELAVANCIN', 'CCCCCCCCCCNCCN[C@@]1(C)C[C@H](O[C@H]2[C@H](Oc3c4cc5cc3Oc3ccc(cc3Cl)[C@@H](O)[C@@H](NC(=O)[C@@H](CC(C)C)NC)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@H]5C(=O)N[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)O)c5cc(O)c(CNCP(=O)(O)O)c(O)c5-c5cc3ccc5O)[C@H](O)c3ccc(c(Cl)c3)O4)O[C@H](CO)[C@@H](O)[C@@H]2O)O[C@@H](C)[C@H]1O',
     'Protein', 'MOL', 'macrocyclic_peptide'),
    # OCTREOTIDE: 
    ('OCTREOTIDE', 'C[C@@H](O)[C@@H]1NC(=O)[C@H](CCCCN)NC(=O)[C@@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H](NC(=O)[C@H](N)Cc2ccccc2)CSSC[C@@H](C(=O)N[C@H](CO)[C@@H](C)O)NC1=O',
     'Protein', 'MOL', 'macrocyclic_peptide'),
    # CARFILZOMIB: tetrapeptide epoxyketone, MW 720, linear -- the shortest true peptide in the set.
    ('CARFILZOMIB', 'CC(C)C[C@H](NC(=O)[C@H](CCc1ccccc1)NC(=O)CN1CCOCC1)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(C)C)C(=O)[C@@]1(C)CO1',
     'Protein', 'BOTH', 'peptide'),
    # BACITRACIN: 
    ('BACITRACIN', 'CC[C@H](C)[C@H](N)C1=N[C@H](C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CCC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CCCCNC(=O)[C@H](CC(N)=O)NC(=O)[C@@H](CC(=O)O)NC(=O)[C@H](Cc3cnc[nH]3)NC(=O)[C@@H](Cc3ccccc3)NC(=O)[C@H]([C@@H](C)CC)NC(=O)[C@@H](CCCN)NC2=O)[C@@H](C)CC)CS1',
     'Protein', 'BOTH', 'macrocyclic_peptide'),
    # LEUPROLIDE: 
    ('LEUPROLIDE', 'CCNC(=O)[C@@H]1CCCN1C(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CO)NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@@H]1CCC(=O)N1',
     'Protein', 'MOL', 'peptide'),
    # VANCOMYCIN: glycopeptide antibiotic; ChEMBL types it `Small molecule`. Backbone fraction 0.248 is below the peptide floor -- caught by the glycopeptide gate.
    ('VANCOMYCIN', 'CN[C@H](CC(C)C)C(=O)N[C@H]1C(=O)N[C@@H](CC(N)=O)C(=O)N[C@H]2C(=O)N[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)O)c4cc(O)cc(O)c4-c4cc3ccc4O)[C@H](O)c3ccc(c(Cl)c3)Oc3cc2cc(c3O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O[C@H]2C[C@](C)(N)[C@H](O)[C@H](C)O2)Oc2ccc(cc2Cl)[C@H]1O',
     'Small molecule', 'MOL', 'macrocyclic_peptide'),
    # CYCLOSPORINE: 
    ('CYCLOSPORINE', 'C/C=C/C[C@@H](C)[C@@H](O)[C@H]1C(=O)N[C@@H](CC)C(=O)N(C)CC(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](CC(C)C)C(=O)N(C)[C@@H](C(C)C)C(=O)N1C',
     'Protein', 'BOTH', 'macrocyclic_peptide'),
    # CHEMBL5394125: IL-17C macrocyclic peptide; ChEMBL molecule_type is NULL for all 98.
    ('CHEMBL5394125', 'CC(C)C[C@@H]1NC(=O)[C@H](CO)NC(=O)CNC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](C(C)C)NC(=O)[C@H](Cc2cnc[nH]2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)CSC[C@@H](C(=O)NCC(N)=O)NC(=O)CNC1=O',
     'NULL', 'MOL', 'macrocyclic_peptide'),
    # CHEMBL5394176: IL-17C macrocyclic peptide; ChEMBL molecule_type is NULL for all 98.
    ('CHEMBL5394176', 'CC[C@H](C)[C@@H]1NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)CSC[C@@H](C(=O)NCC(N)=O)NC(=O)CNC(=O)[C@H](C)NC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CC(C)C)NC1=O',
     'NULL', 'MOL', 'macrocyclic_peptide'),
    # CHEMBL5620040: IL-17C macrocyclic peptide; ChEMBL molecule_type is NULL for all 98.
    ('CHEMBL5620040', 'CC[C@H](C)[C@@H]1NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)CSC[C@@H](C(=O)NCC(N)=O)NC(=O)CNC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CC(C)C)NC1=O',
     'NULL', 'MOL', 'macrocyclic_peptide'),
    # CHEMBL5620164: IL-17C macrocyclic peptide; ChEMBL molecule_type is NULL for all 98.
    ('CHEMBL5620164', 'CC[C@H](C)[C@@H]1NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)CSC[C@@H](C(=O)NCC(N)=O)NC(=O)CNC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CC(C)C)NC1=O',
     'NULL', 'MOL', 'macrocyclic_peptide'),
    # CHEMBL5091495: IL-17A difluorocyclohexyl series -- the potent compounds molecule_type abstains on.
    ('CHEMBL5091495', 'CC(F)(F)CNC(=O)[C@@H](CC(F)F)c1ccc2[nH]c([C@@H](NC(=O)c3nonc3C3CC3)C3CCC(F)(F)CC3)nc2c1F',
     'Unknown', 'MOL', 'small_molecule'),
    # CHEMBL4870813: IL-17A difluorocyclohexyl series -- the potent compounds molecule_type abstains on.
    ('CHEMBL4870813', 'Cc1nonc1C(=O)N[C@H](c1nc2c(F)c([C@H](CC(F)(F)F)c3nnc4ccc(C#N)cn34)ccc2[nH]1)C1CCC(F)(F)CC1',
     'Unknown', 'MOL', 'small_molecule'),
    # CHEMBL5080633: IL-17A difluorocyclohexyl series -- the potent compounds molecule_type abstains on.
    ('CHEMBL5080633', 'CC(F)(F)CNC(=O)C(CC(F)F)c1ccc2[nH]c([C@@H](NC(=O)c3nonc3C3CC3)C3CCC(F)(F)CC3)nc2c1F',
     'Unknown', 'MOL', 'small_molecule'),
    # CHEMBL5081660: IL-17A difluorocyclohexyl series -- the potent compounds molecule_type abstains on.
    ('CHEMBL5081660', 'CNC(=O)[C@H](NC(=O)[C@@H](CC(F)(F)F)c1ccc2[nH]c([C@@H](NC(=O)c3nonc3C)C3CCC(F)(F)CC3)nc2c1F)C(C)C',
     'Unknown', 'MOL', 'small_molecule'),
    # CHEMBL4854787: 
    ('CHEMBL4854787', 'CC(C)C[C@H](NC(=O)[C@H](C)NC(=O)[C@@H](N)CS)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CS)C(=O)O)C(C)C)[C@@H](C)O',
     'Unknown', 'MOL', 'peptide'),
    # CHEMBL5069442: 
    ('CHEMBL5069442', 'CSCC[C@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@@H](N)CS)C(C)C)C(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@@H](CS)C(=O)N[C@@H](CCCNC(=N)N)C(=O)O',
     'Unknown', 'MOL', 'peptide'),
    # CHEMBL4860540: 
    ('CHEMBL4860540', 'CC(C)C[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](N)CSSC[C@@H](C(N)=O)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](C(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccc(O)cc2)NC(=O)CNC1=O',
     'Unknown', 'MOL', 'macrocyclic_peptide'),
    # CHEMBL4869025: 
    ('CHEMBL4869025', 'CC(C)C[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](NC(=O)[C@H](CO)NC(=O)CN)CSSC[C@@H](C(=O)NCC(=O)NCC(=O)NCC(=O)N[C@@H](CO)C(=O)N[C@@H](Cc2cnc[nH]2)C(=O)N[C@@H](Cc2cnc[nH]2)C(=O)N[C@@H](Cc2cnc[nH]2)C(=O)N[C@@H](Cc2cnc[nH]2)C(=O)N[C@@H](Cc2cnc[nH]2)C(=O)N[C@@H](Cc2cnc[nH]2)C(=O)O)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](C(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccc(O)cc2)NC(=O)CNC1=O',
     'Unknown', 'MOL', 'macrocyclic_peptide'),
    # EVERNIMICIN: ChEMBL molecule_type = Oligosaccharide.
    ('EVERNIMICIN', 'COC[C@H]1O[C@@H](O[C@@H]2OC[C@@H]3O[C@@]4(OC[C@@H](OC(=O)c5c(C)cc(O)cc5O)[C@@H]5OCO[C@H]54)O[C@H]3[C@H]2O)[C@@H](OC)[C@@H](O)[C@@H]1O[C@@H]1O[C@H](C)[C@H](OC)[C@H](O[C@@H]2O[C@H](C)[C@H]3O[C@]4(C[C@@H](O)[C@H](O[C@H]5C[C@@H](O[C@H]6C[C@](C)([N+](=O)[O-])[C@@H](OC)[C@H](C)O6)[C@H](OC(=O)c6c(C)c(Cl)c(O)c(Cl)c6OC)[C@@H](C)O5)[C@@H](C)O4)O[C@]3(C)[C@@H]2O)[C@H]1O',
     'Oligosaccharide', 'MOL', 'oligosaccharide'),
    # GLYCYRRHIZIN: ChEMBL molecule_type = Oligosaccharide.
    ('GLYCYRRHIZIN', 'CC1(C)[C@@H](O[C@H]2O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]2O[C@@H]2O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]2O)CC[C@]2(C)[C@H]3C(=O)C=C4[C@@H]5C[C@@](C)(C(=O)O)CC[C@]5(C)CC[C@@]4(C)[C@]3(C)CC[C@@H]12',
     'Oligosaccharide', 'MOL', 'oligosaccharide'),
    # LIVIDOMYCIN: ChEMBL molecule_type = Oligosaccharide.
    ('LIVIDOMYCIN', 'NC[C@@H]1O[C@H](O[C@H]2[C@@H](O)[C@H](O[C@@H]3[C@@H](O)[C@H](N)C[C@H](N)[C@H]3O[C@H]3O[C@H](CO)[C@@H](O)C[C@H]3N)O[C@@H]2CO)[C@H](N)[C@@H](O)[C@@H]1O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O',
     'Oligosaccharide', 'MOL', 'oligosaccharide'),
    # ECHINACOSIDE: ChEMBL molecule_type = Oligosaccharide.
    ('ECHINACOSIDE', 'C[C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](OCCc3ccc(O)c(O)c3)O[C@H](CO[C@H]3O[C@@H](CO)[C@H](O)[C@@H](O)[C@@H]3O)[C@H]2OC(=O)/C=C/c2ccc(O)c(O)c2)[C@@H](O)[C@@H](O)[C@@H]1O',
     'Oligosaccharide', 'MOL', 'oligosaccharide'),
]


PEPTIDE_FAMILY = {"peptide", "macrocyclic_peptide", "protein_or_antibody"}


def _coarse(m):
    return "peptide_family" if m in PEPTIDE_FAMILY else m


def main() -> int:
    from modality import classify_compound

    fails, coarse_fails = [], []
    for cid, smi, mt, st, truth in CASES:
        mt = None if mt == "NULL" else mt
        st = None if st == "NULL" else st
        v = classify_compound(smi, mt, st)
        if v.modality != truth:
            fails.append((cid, truth, v.modality, v.reason))
        if _coarse(v.modality) != _coarse(truth):
            coarse_fails.append(cid)

    n = len(CASES)
    print(f"modality selftest: {n} cases")
    print(f"  exact:  {n - len(fails)}/{n}")
    print(f"  coarse: {n - len(coarse_fails)}/{n}   (small_molecule vs not)")

    dangerous = [f for f in fails if f[2] == "small_molecule"]
    print(f"  FALSE small-molecule calls: {len(dangerous)}")

    for cid, truth, got, why in fails:
        print(f"    MISS {cid}: expected {truth}, got {got}\n         {why}")

    if coarse_fails or dangerous:
        print("FAILED")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
