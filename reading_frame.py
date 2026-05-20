from information_transfer import BioInfoTrans


class ReadingFrameAnalysis:
    def __init__(self, bioinfo: BioInfoTrans):
        self.bioinfo = bioinfo

    def three_frames_translation(self, bases, RNA_polarity=["5'", "3'"]):
        frame_1 = self.bioinfo.match_codons(
            bases, RNA_polarity, name_type="Letter", frame=1, terminus=False
        )
        frame_2 = self.bioinfo.match_codons(
            bases, RNA_polarity, name_type="Letter", frame=2, terminus=False
        )
        frame_3 = self.bioinfo.match_codons(
            bases, RNA_polarity, name_type="Letter", frame=3, terminus=False
        )

        protein_1 = "  ".join(frame_1)
        protein_2 = " " + "  ".join(frame_2)
        protein_3 = "  " + "  ".join(frame_3)

        print(f"{protein_3}\n{protein_2}\n{protein_1}\n{bases}")

    def six_frames_translation(self, bases, RNA_polarity=["5'", "3'"]):
        frame_1 = self.bioinfo.match_codons(
            bases, RNA_polarity, name_type="Letter", frame=1, terminus=False
        )
        frame_2 = self.bioinfo.match_codons(
            bases, RNA_polarity, name_type="Letter", frame=2, terminus=False
        )
        frame_3 = self.bioinfo.match_codons(
            bases, RNA_polarity, name_type="Letter", frame=3, terminus=False
        )

        protein_1 = "  ".join(frame_1)
        protein_2 = " " + "  ".join(frame_2)
        protein_3 = "  " + "  ".join(frame_3)

        c_bases = self.bioinfo.base_pair_RNA(bases)

        c_frame_1 = self.bioinfo.match_codons(
            c_bases, RNA_polarity[::-1], name_type="Letter", frame=1, terminus=False
        )
        c_frame_2 = self.bioinfo.match_codons(
            c_bases, RNA_polarity[::-1], name_type="Letter", frame=3, terminus=False
        )
        c_frame_3 = self.bioinfo.match_codons(
            c_bases, RNA_polarity[::-1], name_type="Letter", frame=2, terminus=False
        )

        c_protein_3 = "   " + "  ".join(c_frame_3)
        c_protein_2 = "    " + "  ".join(c_frame_2)
        c_protein_1 = "  " + "  ".join(c_frame_1)

        print(
            f"{protein_3}\n{protein_2}\n{protein_1}\n{bases}\n{c_bases}\n{c_protein_1}\n{c_protein_2}\n{c_protein_3}"
        )
