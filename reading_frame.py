from information_transfer import BioInfoTrans
from sequence_utils import DNA_or_RNA, parse_polarity


class ReadingFrameAnalysis:
    """
    Analyse protein translation across multiple reading frames.

    Supports analysis of:
    - 3-frame translation (single strand)
    - 6-frame translation (both strands with complementary strand)

    Attributes:
        bioinfo (BioInfoTrans): BioInfoTrans instance for codon translation.
        transcribed_rna (str): RNA strand with polarity markers (e.g., "5'AUGC3'").
        rna_bases (str): RNA sequence without polarity markers.
        rna_polarity (Tuple[str, str]): (left_polarity, right_polarity) e.g., ("5'", "3'").
    """

    def __init__(self, bioinfo: BioInfoTrans, transcribed_rna: str) -> None:
        """
        Initialize ReadingFrameAnalysis with RNA strand.

        Args:
            bioinfo (BioInfoTrans): BioInfoTrans instance for codon lookup.
            transcribed_rna (str): RNA sequence with polarity markers (e.g., "5'AUGC3'").

        Raises:
            ValueError: If strand lacks proper polarity or is not valid RNA.

        Examples:
            >>> analysis = ReadingFrameAnalysis(bio, "5'AUGC3'")
        """
        self.bioinfo = bioinfo
        self.transcribed_rna = transcribed_rna

        self.validate_rna_strand()

    def validate_rna_strand(self) -> None:
        """
        Validate the RNA strand for correct polarity and valid RNA bases.

        Checks:
        - Strand has proper 5'/3' polarity markers at both ends (if none, defaults to 5'-3')
        - Sequence contains only valid RNA bases (AUGC)

        Sets:
        - self.rna_polarity (Tuple[str, str]): The detected polarity markers (left_polarity, right_polarity).
        - self.rna_bases (str): Extracted RNA sequence without markers

        Raises:
            ValueError: If polarity is missing or strand is not valid RNA.
        """
        self.rna_bases, self.rna_polarity = parse_polarity(self.transcribed_rna)

        # check the strand is RNA
        _, RNA = DNA_or_RNA(self.rna_bases)
        if not RNA:
            raise ValueError("Strand is not RNA")

        # for visualization, we want the bases to be displayed 5' on left, 3' on right
        if self.rna_polarity == ("3'", "5'"):
            self.rna_bases = self.rna_bases[::-1]
            self.rna_polarity = self.rna_polarity[::-1]

    def three_frames_translation(self) -> None:
        """
        Display translation of RNA in all three possible reading frames.

        Prints protein sequences from frames 1, 2, and 3 aligned with the RNA bases.
        Frame 1 starts at position 0, frame 2 at position 1, frame 3 at position 2.

        Format:
            Line 1: Frame 3 amino acids
            Line 2: Frame 2 amino acids
            Line 3: Frame 1 amino acids
            Line 4: RNA sequence

        Returns:
            Tuple[List[str], List[str], List[str]]: A tuple containing the three frames as lists of amino acids.

        Examples:
            >>> analysis.three_frames_translation()
            (prints formatted frame analysis)
        """
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

        return frame_1, frame_2, frame_3

    def six_frames_translation(self) -> None:
        """
        Display translation of RNA in all six possible reading frames.

        Analyses both the forward strand (3 frames) and the complementary strand (3 frames).
        Complementary strand uses reverse polarity for proper frame alignment.

        Prints protein sequences in order:
            Lines 1-3: Forward strand frames 3, 2, 1
            Line 4: Forward RNA bases
            Line 5: Complementary RNA bases
            Lines 6-8: Complementary strand frames 1, 2, 3

        Returns:
            Tuple[Tuple[List[str], List[str], List[str]], Tuple[List[str], List[str], List[str]]]: A tuple containing the forward strand frames and complementary strand frames as lists of amino acids.

        Examples:
            >>> analysis.six_frames_translation()
            (prints formatted 6-frame analysis)
        """
        # Forward strand frames (5' to 3' RNA to N to C protein)
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

        protein_1 = "+1  " + "  ".join(frame_1)
        protein_2 = "+2   " + "  ".join(frame_2)
        protein_3 = "+3    " + "  ".join(frame_3)

        # Complementary strand frames (3' to 5' RNA to C to N protein)
        c_transcribed_rna = self.bioinfo.base_pair_RNA(self.rna_bases)

        c_frame_1 = self.bioinfo.match_codons(
            c_transcribed_rna,
            self.rna_polarity[::-1],
            name_type="Letter",
            frame=1,
            terminus=False,
        )[::-1]
        c_frame_2 = self.bioinfo.match_codons(
            c_transcribed_rna,
            self.rna_polarity[::-1],
            name_type="Letter",
            frame=2,
            terminus=False,
        )[::-1]
        c_frame_3 = self.bioinfo.match_codons(
            c_transcribed_rna,
            self.rna_polarity[::-1],
            name_type="Letter",
            frame=3,
            terminus=False,
        )[::-1]

        c_protein_3 = "-3     " + "  ".join(c_frame_3)
        c_protein_2 = "-2      " + "  ".join(c_frame_2)
        c_protein_1 = "-1    " + "  ".join(c_frame_1)

        print(
            f"{protein_3}\n{protein_2}\n{protein_1}\n→   {self.rna_bases}\n←   {c_transcribed_rna}\n{c_protein_1}\n{c_protein_2}\n{c_protein_3}"
        )

        return (frame_1, frame_2, frame_3), (c_frame_1, c_frame_2, c_frame_3)
