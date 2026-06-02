"""
Comprehensive unit tests for biological information transfer module.

Tests cover:
- DNA_or_RNA function: validation, edge cases
- gc_content function: calculation, edge cases
- parse_polarity function: strand parsing with polarity markers
- DNA_to_RNA and RNA_to_DNA functions: sequence conversion
- BioInfoTrans class: DNA/RNA conversion, base pairing, replication, transcription, translation
- ReadingFrameAnalysis class: frame analysis and validation
"""

import pytest
import pandas as pd
import warnings
from sequence_utils import (
    DNA_or_RNA,
    gc_content,
    parse_polarity,
    DNA_to_RNA,
    RNA_to_DNA,
)
from information_transfer import BioInfoTrans
from reading_frame import ReadingFrameAnalysis


@pytest.fixture
def codon_table():
    """Load codon table for testing."""
    return pd.read_csv("codon_tables/codon-table-standard-rna.csv")


class TestDNAorRNA:
    """Tests for DNA_or_RNA function."""

    def test_valid_dna_only(self):
        """DNA sequence should return (True, False)."""
        assert DNA_or_RNA("ATGC") == (True, False)
        assert DNA_or_RNA("TTTT") == (True, False)
        assert DNA_or_RNA("AT") == (True, False)  # Has T, so DNA only

    def test_valid_rna_only(self):
        """RNA sequence should return (False, True)."""
        assert DNA_or_RNA("AUGC") == (False, True)
        assert DNA_or_RNA("UUUU") == (False, True)
        assert DNA_or_RNA("AU") == (False, True)  # Has U, so RNA only

    def test_valid_both_dna_and_rna(self):
        """Sequence with only A, G, C should return (True, True)."""
        assert DNA_or_RNA("AGC") == (True, True)
        assert DNA_or_RNA("AAAGGGCCC") == (True, True)

    def test_invalid_character(self):
        """Sequence with invalid bases should return (False, False)."""
        assert DNA_or_RNA("ATGCX") == (False, False)
        assert DNA_or_RNA("N") == (False, False)

    def test_both_t_and_u(self):
        """Sequence with both T and U should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot contain both T and U"):
            DNA_or_RNA("ATGU")

    def test_empty_sequence(self):
        """Empty sequence should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            DNA_or_RNA("")

    def test_case_insensitivity(self):
        """Function should handle lowercase input."""
        assert DNA_or_RNA("atgc") == (True, False)
        assert DNA_or_RNA("augc") == (False, True)
        assert DNA_or_RNA("AtGc") == (True, False)

    def test_single_base_dna(self):
        """Single DNA base should be recognized."""
        assert DNA_or_RNA("A") == (True, True)  # A is valid in both
        assert DNA_or_RNA("T") == (True, False)
        assert DNA_or_RNA("G") == (True, True)  # G is valid in both
        assert DNA_or_RNA("C") == (True, True)  # C is valid in both

    def test_single_base_rna(self):
        """Single RNA base should be recognized."""
        assert DNA_or_RNA("U") == (False, True)


class TestGCContent:
    """Tests for gc_content function."""

    def test_fifty_percent_gc(self):
        """Sequence with 50% GC content."""
        assert gc_content("ATGC") == 50.0
        assert gc_content("GGCC") == 100.0

    def test_zero_percent_gc(self):
        """Sequence with no G or C."""
        assert gc_content("AAAA") == 0.0
        assert gc_content("TTTT") == 0.0
        assert gc_content("UUUU") == 0.0

    def test_hundred_percent_gc(self):
        """Sequence with 100% GC content."""
        assert gc_content("GCGC") == 100.0
        assert gc_content("GGGG") == 100.0

    def test_single_base_g(self):
        """Single G should be 100% GC."""
        assert gc_content("G") == 100.0

    def test_single_base_c(self):
        """Single C should be 100% GC."""
        assert gc_content("C") == 100.0

    def test_single_base_a(self):
        """Single A should be 0% GC."""
        assert gc_content("A") == 0.0

    def test_single_base_t(self):
        """Single T should be 0% GC."""
        assert gc_content("T") == 0.0

    def test_rna_sequence(self):
        """RNA sequence should calculate GC content correctly."""
        assert gc_content("AUGC") == 50.0
        assert gc_content("UUUU") == 0.0
        assert gc_content("GGGG") == 100.0

    def test_mixed_case(self):
        """Function should handle mixed case."""
        assert gc_content("AtGc") == 50.0
        assert gc_content("AuGc") == 50.0

    def test_long_sequence(self):
        """Test with longer sequence."""
        seq = "ATGCATGCATGC"  # 6 out of 12 are G or C = 50%
        expected = 50.0
        assert abs(gc_content(seq) - expected) < 0.01

    def test_empty_sequence(self):
        """Empty sequence should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            gc_content("")

    def test_invalid_sequence(self):
        """Invalid sequence should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence must be valid DNA or RNA"):
            gc_content("ATGCX")

    def test_both_t_and_u(self):
        """Sequence with both T and U should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence cannot contain both T and U"):
            gc_content("ATGCU")


class TestParsePolarity:
    """Tests for parse_polarity function."""

    def test_valid_5_to_3_polarity(self):
        """Standard 5' to 3' polarity should be parsed correctly."""
        bases, polarity = parse_polarity("5'ATGC3'")
        assert bases == "ATGC"
        assert polarity == ("5'", "3'")

    def test_valid_3_to_5_polarity(self):
        """3' to 5' polarity should be parsed correctly."""
        bases, polarity = parse_polarity("3'ATGC5'")
        assert bases == "ATGC"
        assert polarity == ("3'", "5'")

    def test_no_polarity_provided(self):
        """Strand without polarity should default to 5'-3' with warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bases, polarity = parse_polarity("ATGC")
            assert len(w) == 1
            assert "No strand polarity provided" in str(w[0].message)
            assert bases == "ATGC"
            assert polarity == ("5'", "3'")

    def test_spaces_removed(self):
        """Spaces should be removed during parsing."""
        bases, polarity = parse_polarity("5' A T G C 3'")
        assert bases == "ATGC"
        assert polarity == ("5'", "3'")

    def test_case_insensitive(self):
        """Parsing should handle lowercase bases."""
        bases, polarity = parse_polarity("5'atgc3'")
        assert bases == "ATGC"
        assert polarity == ("5'", "3'")

    def test_invalid_polarity_partial(self):
        """Partial polarity should raise ValueError."""
        with pytest.raises(ValueError, match="Strand does not have correct polarity"):
            parse_polarity("5'ATGC")

    def test_invalid_polarity_wrong_combo(self):
        """Wrong polarity combination should raise ValueError."""
        with pytest.raises(ValueError, match="Strand does not have correct polarity"):
            parse_polarity("5'ATGC5'")

    def test_very_short_sequence(self):
        """Very short sequence should work."""
        bases, polarity = parse_polarity("5'AT3'")
        assert bases == "AT"
        assert polarity == ("5'", "3'")


class TestDNAtoRNAConversion:
    """Tests for DNA_to_RNA function."""

    def test_dna_to_rna_simple(self):
        """DNA sequence should convert T to U."""
        assert DNA_to_RNA("ATGC") == "AUGC"
        assert DNA_to_RNA("TTTT") == "UUUU"

    def test_dna_to_rna_warning_invalid(self):
        """Invalid DNA should trigger warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = DNA_to_RNA("AUGC")  # RNA, not DNA
            assert len(w) == 1
            assert "not a valid DNA sequence" in str(w[0].message)
            assert result == "AUGC"  # Still converts


class TestRNAtoDNAConversion:
    """Tests for RNA_to_DNA function."""

    def test_rna_to_dna_simple(self):
        """RNA sequence should convert U to T."""
        assert RNA_to_DNA("AUGC") == "ATGC"
        assert RNA_to_DNA("UUUU") == "TTTT"

    def test_rna_to_dna_warning_invalid(self):
        """Invalid RNA should trigger warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = RNA_to_DNA("ATGC")  # DNA, not RNA
            assert len(w) == 1
            assert "not a valid RNA sequence" in str(w[0].message)
            assert result == "ATGC"  # Still converts


class TestBioInfoTransInitialization:
    """Tests for BioInfoTrans initialization and validation."""

    def test_valid_initialization_with_polarity(self, codon_table):
        """Valid DNA with polarity markers should initialize."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        assert bio.bases == "ATGC"
        assert bio.template_polarity == ("5'", "3'")

    def test_valid_initialization_without_polarity(self, codon_table):
        """DNA without polarity should initialize with warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            bio = BioInfoTrans("ATGC", codon_table)
            assert len(w) == 1
            assert "No strand polarity provided" in str(w[0].message)
            assert bio.template_polarity == ("5'", "3'")
            assert bio.bases == "ATGC"

    def test_case_insensitive(self, codon_table):
        """Initialization should be case-insensitive."""
        bio1 = BioInfoTrans("5'atgc3'", codon_table)
        bio2 = BioInfoTrans("5'ATGC3'", codon_table)
        assert bio1.bases == bio2.bases

    def test_valid_both_polarities(self, codon_table):
        """Both 5'-3' and 3'-5' polarity orders should be valid."""
        bio1 = BioInfoTrans("5'ATGC3'", codon_table)
        bio2 = BioInfoTrans("3'ATGC5'", codon_table)
        # Both should initialize successfully
        assert bio1.bases == "ATGC"
        assert bio2.bases == "ATGC"

    def test_invalid_polarity_partial(self, codon_table):
        """Partial polarity should raise ValueError."""
        with pytest.raises(ValueError, match="Strand does not have correct polarity"):
            BioInfoTrans("5'ATGC", codon_table)

    def test_non_dna_sequence(self, codon_table):
        """RNA sequence should raise ValueError."""
        with pytest.raises(ValueError, match="Strand is not DNA"):
            BioInfoTrans("5'AUGC3'", codon_table)

    def test_invalid_bases(self, codon_table):
        """Invalid bases should raise ValueError."""
        with pytest.raises(ValueError):
            BioInfoTrans("5'ATGCX3'", codon_table)


class TestBioInfoTransBasePairing:
    """Tests for base pairing methods."""

    def test_base_pair_dna(self, codon_table):
        """DNA base pairing should follow Watson-Crick rules."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        assert bio.base_pair_DNA("ATGC") == "TACG"
        assert bio.base_pair_DNA("AAAA") == "TTTT"
        assert bio.base_pair_DNA("GGGG") == "CCCC"

    def test_base_pair_rna(self, codon_table):
        """RNA base pairing should follow Watson-Crick rules."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        assert bio.base_pair_RNA("AUGC") == "UACG"
        assert bio.base_pair_RNA("AAAA") == "UUUU"
        assert bio.base_pair_RNA("GGGG") == "CCCC"


class TestBioInfoTransMatchCodons:
    """Tests for codon matching in translation from RNA to protein."""

    def test_match_codons_valid(self, codon_table):
        """Valid codon should return correct amino acid."""
        bio = BioInfoTrans("3' TACTTACACACG 5'", codon_table)
        protein = bio.match_codons(
            "AUGAAUGUGUGC", ("5'", "3'"), name_type="Short Name", frame=1
        )
        assert protein == ["N", "Met", "Asn", "Val", "Cys", "C"]

    def test_match_codons_invalid(self, codon_table):
        """Invalid codon should raise ValueError."""
        bio = BioInfoTrans("3' TACTTACACACG 5'", codon_table)
        with pytest.raises(ValueError, match="Unknown codon: XXX"):
            bio.match_codons("AUGXXX", ("5'", "3'"), name_type="Short Name", frame=1)

    def test_match_codons_polarity_mirror(self, codon_table):
        """Mirrored polarity should be handled correctly - protein sequence should be the same."""
        bio = BioInfoTrans("3' TACTTACACACG 5'", codon_table)
        protein1 = bio.match_codons(
            "AUGAAU", ("5'", "3'"), name_type="Short Name", frame=1
        )
        protein2 = bio.match_codons(
            "UAAGUA", ("3'", "5'"), name_type="Short Name", frame=1
        )
        assert protein1 == ["N", "Met", "Asn", "C"]
        assert protein2 == ["N", "Met", "Asn", "C"]

    def test_match_codons_non_multiple_of_three(self, codon_table):
        """Non-multiple-of-three RNA sequence should ignore incomplete trailing codons."""
        bio = BioInfoTrans("3' TACTTACACACG 5'", codon_table)
        protein = bio.match_codons(
            "AUGAAUGUGUG", ("5'", "3'"), name_type="Short Name", frame=1
        )
        assert protein == [
            "N",
            "Met",
            "Asn",
            "Val",
            "C",
        ]  # Last 'UG' is incomplete and should be ignored

    def test_match_codons_frame_shift(self, codon_table):
        """Frame shift should return different amino acid."""
        bio = BioInfoTrans("3' TACTTACACACG 5'", codon_table)
        protein_frame1 = bio.match_codons(
            "AUGUCUACUGCCGGAUGGCUAAUAUAA", ("5'", "3'"), name_type="Short Name", frame=1
        )
        protein_frame2 = bio.match_codons(
            "AUGUCUACUGCCGGAUGGCUAAUAUAA",
            ("5'", "3'"),
            name_type="Short Name",
            frame=2,
        )
        protein_frame3 = bio.match_codons(
            "AUGUCUACUGCCGGAUGGCUAAUAUAA", ("5'", "3'"), name_type="Short Name", frame=3
        )
        assert protein_frame1 == [
            "N",
            "Met",
            "Ser",
            "Thr",
            "Ala",
            "Gly",
            "Trp",
            "Leu",
            "Ile",
            "Stop",
            "C",
        ]
        assert protein_frame2 == [
            "N",
            "Cys",
            "Leu",
            "Leu",
            "Pro",
            "Asp",
            "Gly",
            "Stop",
            "Tyr",
            "C",
        ]
        assert protein_frame3 == [
            "N",
            "Val",
            "Tyr",
            "Cys",
            "Arg",
            "Met",
            "Ala",
            "Asn",
            "Ile",
            "C",
        ]


class TestBioInfoTransReplication:
    """Tests for DNA replication."""

    def test_replicate_dna_simple(self, codon_table):
        """DNA replication should produce complementary strand."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        replicated = bio.replicate_DNA()
        assert replicated == "3'TACG5'"

    def test_replicate_dna_polarity_swap(self, codon_table):
        """Replication should swap polarity."""
        bio = BioInfoTrans("5'AAA3'", codon_table)
        replicated = bio.replicate_DNA()
        assert replicated.startswith("3'")
        assert replicated.endswith("5'")


class TestBioInfoTransTranscription:
    """Tests for DNA to RNA transcription."""

    def test_transcribe_simple(self, codon_table):
        """Transcription should produce complementary RNA."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        transcribed = bio.transcribe_RNA()
        assert transcribed == "3'UACG5'"

    def test_transcribe_contains_uracil(self, codon_table):
        """Transcribed RNA should contain U, not T."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        transcribed = bio.transcribe_RNA()
        assert "U" in transcribed
        assert "T" not in transcribed


class TestBioInfoTransTranslation:
    """Tests for RNA to protein translation."""

    def test_translate_simple_start_codon(self, codon_table):
        """Translation should start with methionine (M)."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        proteins = bio.translate_RNA(
            "AUGC", ("5'", "3'"), name_type="Letter", terminus=False
        )
        assert proteins[0] == "M"  # AUG codes for Met (start)

    def test_translate_with_terminus(self, codon_table):
        """Translation with terminus should add N and C markers."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        proteins = bio.translate_RNA(
            "AUGC", ("5'", "3'"), name_type="Letter", terminus=True
        )
        assert proteins[0] == "N"  # N terminus at start
        assert proteins[-1] == "C"  # C terminus at end

    def test_translate_without_terminus(self, codon_table):
        """Translation without terminus should not add N/C."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        proteins = bio.translate_RNA(
            "AUGC", ("5'", "3'"), name_type="Letter", terminus=False
        )
        assert proteins[0] == "M"
        assert proteins[-1] != "C"


class TestBioInfoTransInfoTransfer:
    """Tests for complete information transfer workflow."""

    def test_info_transfer_returns_tuple(self, codon_table, capsys):
        """info_transfer should return tuple of (DNA, RNA, protein)."""
        bio = BioInfoTrans("5'ATGATGATG3'", codon_table)
        result = bio.info_transfer()
        assert isinstance(result, tuple)
        assert len(result) == 3

    def test_info_transfer_prints_output(self, codon_table, capsys):
        """info_transfer should print intermediate steps."""
        bio = BioInfoTrans("5'ATGATGATG3'", codon_table)
        bio.info_transfer()
        captured = capsys.readouterr()
        assert "DNA Template Strand" in captured.out
        assert "DNA Replication" in captured.out
        assert "RNA Transcription" in captured.out
        assert "Protein Translation" in captured.out


class TestReadingFrameAnalysisInitialization:
    """Tests for ReadingFrameAnalysis initialization."""

    def test_valid_initialization(self, codon_table):
        """Valid RNA strand should initialize."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        transcribed = bio.transcribe_RNA()
        analysis = ReadingFrameAnalysis(bio, transcribed)
        assert analysis.bioinfo == bio

    def test_missing_polarity(self, codon_table):
        """RNA without polarity should initialize with warning and should default to 5'-3'."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            analysis = ReadingFrameAnalysis(bio, "AUGC")
            assert len(w) == 1
            assert "No strand polarity provided" in str(w[0].message)
            assert analysis.rna_polarity == ("5'", "3'")
            assert analysis.rna_bases == "AUGC"

    def test_invalid_polarity(self, codon_table):
        """RNA with invalid polarity should raise ValueError."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        with pytest.raises(ValueError, match="Strand does not have correct polarity"):
            ReadingFrameAnalysis(bio, "5'AUGC")

    def test_dna_instead_of_rna(self, codon_table):
        """DNA sequence should raise ValueError."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        with pytest.raises(ValueError, match="Strand is not RNA"):
            ReadingFrameAnalysis(bio, "5'ATGC3'")


class TestReadingFrameAnalysisValidation:
    """Tests for ReadingFrameAnalysis validation."""

    def test_rna_bases_extracted_correctly(self, codon_table):
        """RNA bases should be extracted without polarity markers."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        transcribed = bio.transcribe_RNA()
        analysis = ReadingFrameAnalysis(bio, transcribed)
        assert "5'" not in analysis.rna_bases
        assert "3'" not in analysis.rna_bases
        assert len(analysis.rna_bases) > 0

    def test_rna_polarity_stored_correctly(self, codon_table):
        """RNA polarity should be stored as tuple."""
        bio = BioInfoTrans("5'ATGC3'", codon_table)
        transcribed = bio.transcribe_RNA()
        analysis = ReadingFrameAnalysis(bio, transcribed)
        assert isinstance(analysis.rna_polarity, tuple)
        assert len(analysis.rna_polarity) == 2


class TestReadingFrameAnalysisTranslation:
    """Tests for reading frame translation methods."""

    def test_three_frames_translation_prints(self, codon_table, capsys):
        """three_frames_translation should print output."""
        bio = BioInfoTrans("5'ATGATGATG3'", codon_table)
        transcribed = bio.transcribe_RNA()
        analysis = ReadingFrameAnalysis(bio, transcribed)
        analysis.three_frames_translation()
        captured = capsys.readouterr()
        assert len(captured.out) > 0

    def test_three_frames_translation_returns_tuple(self, codon_table):
        """three_frames_translation should return tuple of frames."""
        bio = BioInfoTrans("5'ATGATGATG3'", codon_table)
        analysis = ReadingFrameAnalysis(bio, "5'AUGAAUGUGUGC3'")
        forward_frames = analysis.three_frames_translation()
        assert forward_frames == (
            ["M", "N", "V", "C"],
            ["*", "M", "C"],
            ["E", "C", "V"],
        )

    def test_six_frames_translation_prints(self, codon_table, capsys):
        """six_frames_translation should print output."""
        bio = BioInfoTrans("5'ATGATGATG3'", codon_table)
        transcribed = bio.transcribe_RNA()
        analysis = ReadingFrameAnalysis(bio, transcribed)
        analysis.six_frames_translation()
        captured = capsys.readouterr()
        assert len(captured.out) > 0

    def test_six_frames_translation_returns_tuple(self, codon_table):
        """six_frames_translation should return tuple of frames."""
        bio = BioInfoTrans("5'ATGATGATG3'", codon_table)
        analysis = ReadingFrameAnalysis(bio, "5'AUGAAUGUGUGC3'")
        forward_frames, reverse_frames = analysis.six_frames_translation()
        assert forward_frames == (
            ["M", "N", "V", "C"],
            ["*", "M", "C"],
            ["E", "C", "V"],
        )
        assert reverse_frames == (
            ["H", "I", "H", "A"],
            ["F", "T", "H"],
            ["S", "H", "T"],
        )


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_consecutive_gc_content(self):
        """Consecutive G and C should calculate correctly."""
        assert gc_content("GC" * 50) == 100.0

    def test_alternating_gc_at_content(self):
        """Alternating GC and AT should calculate correctly."""
        assert gc_content("GCAT" * 25) == 50.0

    def test_large_sequence_dna_or_rna(self):
        """Large sequence validation should work."""
        large_dna = "ATGC" * 1000
        is_dna, is_rna = DNA_or_RNA(large_dna)
        assert is_dna
        assert not is_rna

    def test_long_dna_sequence(self, codon_table):
        """Long DNA sequence should initialize."""
        long_seq = "5'" + "ATGC" * 100 + "3'"
        bio = BioInfoTrans(long_seq, codon_table)
        assert len(bio.bases) == 400


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
