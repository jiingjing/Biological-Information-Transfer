from typing import Tuple
import warnings


def DNA_or_RNA(bases: str) -> Tuple[bool, bool]:
    """
    Determine if a nucleotide sequence is valid DNA, RNA, or both.

    A valid DNA sequence contains only ATGC.
    A valid RNA sequence contains only AUGC.

    Args:
        bases (str): The nucleotide sequence to validate.

    Returns:
        Tuple[bool, bool]: A tuple (is_DNA, is_RNA) where:
            - is_DNA is True if the sequence is valid DNA
            - is_RNA is True if the sequence is valid RNA

    Raises:
        ValueError: If sequence is empty or contains both T and U.

    Examples:
        >>> DNA_or_RNA("ATGC")
        (True, False)
        >>> DNA_or_RNA("AUGC")
        (False, True)
        >>> DNA_or_RNA("ATGCU")
        Raises ValueError
    """
    if not bases:
        raise ValueError("Sequence cannot be empty.")

    bases = bases.upper()  # convert to uppercase for case insensitivity

    # If the sequence has both T and U, it's invalid.
    if "T" in bases and "U" in bases:
        raise ValueError("Sequence cannot contain both T and U.")

    DNA_bases = set("ATGC")
    RNA_bases = set("AUGC")

    valid_DNA = set(bases) <= DNA_bases
    valid_RNA = set(bases) <= RNA_bases
    return valid_DNA, valid_RNA


def gc_content(bases: str) -> float:
    """
    Calculate the GC content (percentage of G and C bases) in a nucleotide sequence.

    Args:
        bases (str): A valid DNA or RNA nucleotide sequence.

    Returns:
        float: GC content as a percentage (0-100).

    Raises:
        ValueError: If the sequence is not valid DNA or RNA, or is empty.

    Examples:
        >>> gc_content("ATGC")
        50.0
        >>> gc_content("AAAA")
        0.0
        >>> gc_content("GCGC")
        100.0
    """
    # Validate sequence before processing (also handles empty check)
    is_dna, is_rna = DNA_or_RNA(bases)
    if not (is_dna or is_rna):
        raise ValueError("Sequence must be valid DNA or RNA.")

    bases = bases.upper()  # convert to uppercase for case insensitivity
    GC_count = bases.count("G") + bases.count("C")
    total_count = len(bases)

    return (GC_count / total_count) * 100


def parse_polarity(strand: str) -> Tuple[str, Tuple[str, str]]:
    """
    Parse a nucleotide strand with polarity markers and extract the sequence and polarity.

    Args:
        strand (str): A nucleotide strand with polarity markers (e.g., "5'AUGC3'").

    Returns:
        Tuple[str, Tuple[str, str]]: A tuple containing:
            - The extracted nucleotide sequence without polarity markers, just the bases (e.g., "AUGC")
            - The left polarity marker (e.g., "5'")
            - The right polarity marker (e.g., "3'")

    Raises:
        ValueError: If the strand does not have proper polarity markers.

    Examples:
        >>> parse_polarity("5' AUGC 3'")
        ("AUGC", ("5'", "3'"))
        >>> parse_polarity("3'ATGC5'")
        ("ATGC", ("3'", "5'"))
        >>> parse_polarity("AUGC")
        ("AUGC", ("5'", "3'"))
        >>> parse_polarity("5'AUGC")
        Raises ValueError
    """
    strand = strand.upper().replace(" ", "")  # Remove spaces, convert to uppercase

    left_end = strand[:2]
    right_end = strand[-2:]

    valid_markers = {"3'", "5'"}
    valid_combos = {("5'", "3'"), ("3'", "5'")}

    # Case 1: explicit polarity provided
    if (left_end, right_end) in valid_combos:
        bases = strand[2:-2]
        polarity = (left_end, right_end)
        return bases, polarity

    # Case 2: no polarity provided
    if left_end not in valid_markers and right_end not in valid_markers:
        warnings.warn(
            "No strand polarity provided. Assuming 5'->3'.",
            UserWarning,
            stacklevel=2,
        )
        return strand, ("5'", "3'")

    # Case 3: partial/invalid polarity
    raise ValueError("Strand does not have correct polarity.")


def DNA_to_RNA(bases: str) -> str:
    """
    Convert DNA sequence to RNA by replacing thymine (T) with uracil (U).

    Args:
        bases (str): DNA sequence containing T nucleotides.

    Returns:
        str: RNA sequence with U in place of T.

    Examples:
        >>> DNA_to_RNA("ATGC")
        "AUGC"
    """
    # user warning if strand is not DNA
    DNA, _ = DNA_or_RNA(bases)
    if not DNA:
        warnings.warn("Input is not a valid DNA sequence", UserWarning)
    return bases.replace("T", "U")


def RNA_to_DNA(bases: str) -> str:
    """
    Convert RNA sequence to DNA by replacing uracil (U) with thymine (T).

    Args:
        bases (str): RNA sequence containing U nucleotides.

    Returns:
        str: DNA sequence with T in place of U.

    Examples:
        >>> RNA_to_DNA("AUGC")
        "ATGC"
    """
    # user warning if strand is not RNA
    _, RNA = DNA_or_RNA(bases)
    if not RNA:
        warnings.warn("Input is not a valid RNA sequence", UserWarning)
    return bases.replace("U", "T")
