from information_transfer import BioInfoTrans


class ReadingFrameAnalysis:
    def __init__(self, bioinfo: BioInfoTrans, transcribed_rna):
        self.bioinfo = bioinfo
        self.transcribed_rna = transcribed_rna

        self.validate_rna_strand()

    def validate_rna_strand(self):
        # must be RNA, must have polarity

        ## check polarity - should have '5 or 3' at ends
        left_end = self.transcribed_rna[:2]
        right_end = self.transcribed_rna[-2:]

        if not set(["3'", "5'"]) == set([left_end, right_end]):
            raise ValueError("Strand does not have correct polarity.")

        self.rna_polarity = [left_end, right_end]  # save just the polarity
        self.rna_bases = self.transcribed_rna[2:-2]  # save just the bases

        ## check the strand is RNA
        _, RNA = self.bioinfo.DNA_or_RNA(self.rna_bases)
        if not RNA:
            raise ValueError("Strand is not RNA")

    def three_frames_translation(self):
        frame_1 = self.bioinfo.match_codons(
            self.rna_bases,
            self.rna_polarity,
            name_type="Letter",
            frame=1,
            terminus=False,
        )
        frame_2 = self.bioinfo.match_codons(
            self.rna_bases,
            self.rna_polarity,
            name_type="Letter",
            frame=2,
            terminus=False,
        )
        frame_3 = self.bioinfo.match_codons(
            self.rna_bases,
            self.rna_polarity,
            name_type="Letter",
            frame=3,
            terminus=False,
        )

        protein_1 = "  ".join(frame_1)
        protein_2 = " " + "  ".join(frame_2)
        protein_3 = "  " + "  ".join(frame_3)

        print(f"{protein_3}\n{protein_2}\n{protein_1}\n{self.rna_bases}")

    def six_frames_translation(self):
        frame_1 = self.bioinfo.match_codons(
            self.rna_bases,
            self.rna_polarity,
            name_type="Letter",
            frame=1,
            terminus=False,
        )
        frame_2 = self.bioinfo.match_codons(
            self.rna_bases,
            self.rna_polarity,
            name_type="Letter",
            frame=2,
            terminus=False,
        )
        frame_3 = self.bioinfo.match_codons(
            self.rna_bases,
            self.rna_polarity,
            name_type="Letter",
            frame=3,
            terminus=False,
        )

        protein_1 = "  ".join(frame_1)
        protein_2 = " " + "  ".join(frame_2)
        protein_3 = "  " + "  ".join(frame_3)

        c_transcribed_rna = self.bioinfo.base_pair_RNA(self.rna_bases)

        c_frame_1 = self.bioinfo.match_codons(
            c_transcribed_rna,
            self.rna_polarity[::-1],
            name_type="Letter",
            frame=1,
            terminus=False,
        )
        c_frame_2 = self.bioinfo.match_codons(
            c_transcribed_rna,
            self.rna_polarity[::-1],
            name_type="Letter",
            frame=3,
            terminus=False,
        )
        c_frame_3 = self.bioinfo.match_codons(
            c_transcribed_rna,
            self.rna_polarity[::-1],
            name_type="Letter",
            frame=2,
            terminus=False,
        )

        c_protein_3 = "   " + "  ".join(c_frame_3)
        c_protein_2 = "    " + "  ".join(c_frame_2)
        c_protein_1 = "  " + "  ".join(c_frame_1)

        print(
            f"{protein_3}\n{protein_2}\n{protein_1}\n{self.rna_bases}\n{c_transcribed_rna}\n{c_protein_1}\n{c_protein_2}\n{c_protein_3}"
        )
