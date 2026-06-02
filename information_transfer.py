import warnings
from typing import List, Optional, Tuple
import pandas as pd
from sequence_utils import DNA_or_RNA, DNA_to_RNA, parse_polarity


class BioInfoTrans:
    """
    Central class for biological information transfer between DNA, RNA, and proteins.

    Handles DNA replication, transcription to RNA, and translation to proteins.
    Supports bidirectional strand processing (5'-3' and 3'-5' polarities).

    Attributes:
        template_strand (str): DNA template strand with polarity markers (e.g., "5'ATGC3'").
        bases (str): DNA sequence without polarity markers.
        template_polarity (Tuple[str, str]): [left_polarity, right_polarity] e.g., ["5'", "3'"].
        codon_df (pd.DataFrame): Codon translation lookup table.
    """

    def __init__(self, template_strand: str, codon_df: pd.DataFrame) -> None:
        """
        Initialize BioInfoTrans with a DNA template strand and codon table.

        Args:
            template_strand (str): DNA sequence with polarity markers (e.g., "5'ATGC3'").
                                   Spaces are automatically removed, case-insensitive.
            codon_df (pd.DataFrame): DataFrame with columns "Codon", "Short Name", "Letter"
                                     for codon translation lookup.

        Examples:
            >>> bio = BioInfoTrans("5'ATGC3'", codon_df)
        """
        self.template_strand = template_strand
        self.codon_df = codon_df

        self.validate_template_strand()

    def validate_template_strand(self) -> None:
        """
        Validate the DNA template strand for correct polarity and valid DNA bases.

        Checks:
        - Strand has proper 5'/3' polarity markers at both ends (if none, defaults to 5'-3')
        - Sequence contains only valid DNA bases (ATGC)

        Sets:
        - self.template_polarity (Tuple[str, str]): The detected polarity markers (left_polarity, right_polarity).
        - self.bases (str): Extracted DNA sequence without markers

        Raises:
            ValueError: If polarity is incorrect/incomplete or strand is not valid DNA.

        Warns:
            UserWarning: If no polarity markers are provided (defaults to 5'-3').
        """
        # Parse polarity and extract bases first
        self.bases, self.template_polarity = parse_polarity(self.template_strand)

        # Then validate that bases are DNA
        DNA, _ = DNA_or_RNA(self.bases)
        if not DNA:
            raise ValueError("Strand is not DNA")

    def base_pair_DNA(self, bases: str) -> str:
        """
        Generate complementary DNA strand using Watson-Crick base pairing.

        Pairing rules:
        - A ↔ T
        - G ↔ C

        Args:
            bases (str): DNA sequence.

        Returns:
            str: Complementary DNA strand.

        Examples:
            >>> bio.base_pair_DNA("ATGC")
            "TACG"
        """
        base_table = str.maketrans("ATGC", "TACG")
        base_pair = bases.translate(base_table)
        return base_pair

    def base_pair_RNA(self, bases: str) -> str:
        """
        Generate complementary RNA strand using Watson-Crick base pairing.

        Pairing rules:
        - A ↔ U
        - G ↔ C

        Args:
            bases (str): RNA sequence.

        Returns:
            str: Complementary RNA strand.

        Examples:
            >>> bio.base_pair_RNA("AUGC")
            "UACG"
        """
        base_table = str.maketrans("AUGC", "UACG")
        base_pair = bases.translate(base_table)
        return base_pair

    def match_codons(
        self,
        bases: str,
        rna_polarity: List[str],
        name_type: str = "Short Name",
        frame: int = 1,
        terminus: bool = True,
    ) -> List[str]:
        """
        Translate RNA codons to amino acids based on reading frame and polarity.
        Always returns protein sequence in N to C direction regardless of RNA polarity.
        Ignores incomplete trailing codons.

        Args:
            bases (str): RNA sequence to translate.
            rna_polarity (List[str]): [left_polarity, right_polarity] e.g., ["5'", "3'"].
            name_type (str): Column name in codon_df for output ("Short Name" or "Letter").
                           Defaults to "Short Name".
            frame (int): Reading frame (1, 2, or 3). Frame 1 starts at first base,
                        frame 2 at second, frame 3 at third. Defaults to 1.
            terminus (bool): If True, add N and C terminus markers. Defaults to True.

        Returns:
            List[str]: Translated amino acids or abbreviations in N to C direction.

        Raises:
            ValueError: If codon not found in lookup table.

        Examples:
            >>> proteins = bio.match_codons("AUGCUGA", ("5'", "3'"), frame=1, name_type="Letter", terminus=False)
            ["M", "L"]
        """
        if frame not in (1, 2, 3):
            raise ValueError("Frame must be 1, 2, or 3")

        if rna_polarity[1] == "5'":  # if 3'-5', read right to left i.e. backwards
            bases = bases[
                ::-1
            ]  # reverse the entire RNA sequence for correct codon reading

        start = frame - 1
        triplets = [bases[i : i + 3] for i in range(start, len(bases) - 2, 3)]

        protein = []
        for t in triplets:
            codon = self.codon_df.loc[self.codon_df["Codon"] == t, name_type]
            if len(codon) != 1:
                raise ValueError(f"Unknown codon: {t}")
            protein.append(codon.item())

        # check rna_polarity and insert N/C terminus
        if terminus:
            protein.append("C")  # add C at end
            protein.insert(0, "N")  # add N at start
        return protein

    def replicate_DNA(self) -> str:
        """
        Replicate DNA by generating complementary strand with reversed polarity.

        Returns:
            str: Replicated DNA strand with swapped polarity markers.

        Examples:
            >>> replicated = bio.replicate_DNA()
            "3'TACG5'"
        """
        replicated = self.base_pair_DNA(self.bases)
        return self.template_polarity[1] + replicated + self.template_polarity[0]

    def transcribe_RNA(self) -> str:
        """
        Transcribe DNA to complementary RNA strand.

        Process:
        1. Convert DNA to RNA (T → U)
        2. Generate complementary strand (A↔U, G↔C)
        3. Reverse polarity

        Returns:
            str: Transcribed RNA strand with polarity markers.

        Examples:
            >>> transcribed = bio.transcribe_RNA()
            "3'UACG5'"
        """
        transcribed = DNA_to_RNA(self.bases)
        transcribed = self.base_pair_RNA(transcribed)
        return self.template_polarity[1] + transcribed + self.template_polarity[0]

    def translate_RNA(
        self,
        rna_bases: str,
        rna_polarity: Optional[List[str]] = None,
        name_type: str = "Short Name",
        terminus: bool = True,
    ) -> List[str]:
        """
        Translate RNA sequence to protein amino acids.

        Args:
            rna_bases (str): RNA sequence to translate (without polarity markers).
            rna_polarity (Optional[List[str]]): [left_polarity, right_polarity].
                                              If None, uses reversed template_polarity.
            name_type (str): Column name in codon_df ("Short Name" or "Letter").
            terminus (bool): If True, add N and C terminus markers.

        Returns:
            List[str]: Translated amino acids.

        Examples:
            >>> proteins = bio.translate_RNA("AUGC")
        """
        if rna_polarity is None:
            rna_polarity = self.template_polarity[::-1]
        return self.match_codons(
            rna_bases, rna_polarity, name_type=name_type, terminus=terminus
        )

    def info_transfer(self) -> Tuple[str, str, str]:
        """
        Perform complete central dogma sequence: DNA → RNA → Protein.

        Prints intermediate steps and returns results.

        Returns:
            Tuple[str, str, str]: (replicated_DNA, transcribed_RNA, protein_sequence)

        Examples:
            >>> dna, rna, protein = bio.info_transfer()
        """
        # Template strand
        print(
            f"DNA Template Strand: {self.template_polarity[0] + self.bases + self.template_polarity[1]}"
        )

        # DNA replication
        replicated = self.replicate_DNA()
        print(f"DNA Replication: {replicated}")

        # Transcription
        transcribed = self.transcribe_RNA()
        print(f"RNA Transcription: {transcribed}")

        # Translation (always in N to C direction)
        proteins = self.translate_RNA(rna_bases=transcribed[2:-2])
        proteins = "-".join(proteins)
        print(f"Protein Translation: {proteins}")

        return replicated, transcribed, proteins
