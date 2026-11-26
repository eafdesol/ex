"""Molecule library for the X-H overtone viewer.

Each molecule has one or more forms (anomers / ring tautomers) with an
equilibrium population in water at ~25 C. A form is given either as a
verified isomeric SMILES ("smiles") plus a mapped template, or as a mapped
template plus the CIP descriptors of its stereocentres ("cip", taken from the
IUPAC name), from which the matching stereoisomer is picked automatically.

Atom maps on carbons give labels: 1..9 -> C1..C9, 11..19 -> C1'..C9'.
Populations are rounded literature (NMR) values recalled from memory and not
checked against the original papers.

To add a molecule: add an entry here and run
    python run_xh.py <key>
    python build_xh_viewer.py
"""
from sugar_nir.localmode import SMILES as VERIFIED

T_ALDOPYR = "O[CH:1]1O[CH:5]([CH2:6]O)[CH:4](O)[CH:3](O)[CH:2]1O"
T_MALTOSE = ("O[CH2:16][CH:15]1O[CH:11](O[CH:4]2[CH:5]([CH2:6]O)O[CH:1](O)[CH:2](O)[CH:3]2O)"
             "[CH:12](O)[CH:13](O)[CH:14]1O")
T_KETOPYR = "O[CH2:1][C:2]1(O)O[CH2:6][CH:5](O)[CH:4](O)[CH:3]1O"
T_KETOFUR = "O[CH2:1][C:2]1(O)O[CH:5]([CH2:6]O)[CH:4](O)[CH:3]1O"
T_KETO = "O[CH2:1][C:2](=O)[CH:3](O)[CH:4](O)[CH:5](O)[CH2:6]O"

MOLECULES = {
    "glucose": {
        "name": "포도당", "en": "D-glucose", "kind": "알도헥소스 (단당류)",
        "forms": {
            "α-피라노스": {"smiles": VERIFIED["a-glucose"], "template": T_ALDOPYR, "pop": 0.36},
            "β-피라노스": {"smiles": VERIFIED["b-glucose"], "template": T_ALDOPYR, "pop": 0.64},
        },
    },
    "maltose": {
        "name": "엿당 (말토오스)", "en": "maltose", "kind": "α-1,4 포도당 이당류",
        "n_embed": 400,
        "forms": {
            "α (환원 말단)": {"smiles": VERIFIED["a-maltose"], "template": T_MALTOSE, "pop": 0.42},
            "β (환원 말단)": {"smiles": VERIFIED["b-maltose"], "template": T_MALTOSE, "pop": 0.58},
        },
    },
    "fructose": {
        "name": "과당", "en": "D-fructose", "kind": "케토헥소스 (단당류)",
        "forms": {
            "β-피라노스": {"template": T_KETOPYR, "cip": {2: "R", 3: "S", 4: "R", 5: "R"}, "pop": 0.684},
            "α-피라노스": {"template": T_KETOPYR, "cip": {2: "S", 3: "S", 4: "R", 5: "R"}, "pop": 0.027},
            "β-퓨라노스": {"template": T_KETOFUR, "cip": {2: "R", 3: "S", 4: "S", 5: "R"}, "pop": 0.224},
            "α-퓨라노스": {"template": T_KETOFUR, "cip": {2: "S", 3: "S", 4: "S", 5: "R"}, "pop": 0.062},
            "케토 (열린 사슬)": {"template": T_KETO, "cip": {3: "S", 4: "R", 5: "R"}, "pop": 0.005,
                             "n_embed": 150},
        },
    },
    "allulose": {
        "name": "알룰로스 (D-사이코스)", "en": "D-allulose (D-psicose)", "kind": "케토헥소스, 과당의 C3 에피머",
        "forms": {
            "β-피라노스": {"template": T_KETOPYR, "cip": {2: "R", 3: "R", 4: "R", 5: "R"}, "pop": 0.24},
            "α-피라노스": {"template": T_KETOPYR, "cip": {2: "S", 3: "R", 4: "R", 5: "R"}, "pop": 0.22},
            "β-퓨라노스": {"template": T_KETOFUR, "cip": {2: "R", 3: "R", 4: "S", 5: "R"}, "pop": 0.15},
            "α-퓨라노스": {"template": T_KETOFUR, "cip": {2: "S", 3: "R", 4: "S", 5: "R"}, "pop": 0.39},
            "케토 (열린 사슬)": {"template": T_KETO, "cip": {3: "R", 4: "R", 5: "R"}, "pop": 0.003,
                             "n_embed": 150},
        },
    },
}
