"""Figure numbering for the npj Quantum Information version.

Build scripts and figure files keep their historical names (Fig1-Fig10, FigS1,
FigS2); text is written with those numbers and mapped here to the submitted
numbering. Supplementary figures are numbered in order of appearance in the SI.
"""
import re

MAIN = {1: "Fig. 1", 2: "Fig. 2", 3: "Supplementary Fig. 4", 4: "Fig. 3",
        5: "Supplementary Fig. 3", 6: "Supplementary Fig. 2", 7: "Supplementary Fig. 5",
        8: "Supplementary Fig. 6", 9: "Fig. 4", 10: "Fig. 5"}
SUPP = {1: "Supplementary Fig. 1", 2: "Supplementary Fig. 7"}
FILES = {"Fig1": "Fig1", "Fig2": "Fig2", "Fig4": "Fig3", "Fig9": "Fig4", "Fig10": "Fig5",
         "FigS1": "SuppFig1", "Fig6": "SuppFig2", "Fig5": "SuppFig3", "Fig3": "SuppFig4",
         "Fig7": "SuppFig5", "Fig8": "SuppFig6", "FigS2": "SuppFig7"}
_PAT = re.compile(r"(Supplementary )?Fig\. (\d+)")


def remap(text):
    if not isinstance(text, str):
        return text

    def rep(m):
        n = int(m.group(2))
        return (SUPP if m.group(1) else MAIN)[n]
    return _PAT.sub(rep, text)


def remap_sheet(name):
    """Sheet names such as 'Fig4b_Cisplatin...' or 'Fig10_FigS2_Data'."""
    def rep(m):
        if m.group(1):
            return SUPP[int(m.group(2))].replace("Supplementary Fig. ", "SuppFig")
        return MAIN[int(m.group(2))].replace("Supplementary Fig. ", "SuppFig").replace("Fig. ", "Fig")
    return re.sub(r"Fig(S)?(\d+)", rep, name)
