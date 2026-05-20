class BioInfoTrans:
    def __init__(self, template_strand, codon_df):
        self.template_strand = template_strand  # str: DNA template strand
        self.codon_df = codon_df  # dataframe: codon table

        self.validate_template_strand()

    def validate_template_strand(self):
        # must be DNA, must have polarity

        ## remove any spaces
        template_strand = self.template_strand.replace(" ", "")

        ## check polarity - should have '5 or 3' at ends
        left_end = template_strand[:2]
        right_end = template_strand[-2:]

        if not set(["3'", "5'"]) == set([left_end, right_end]):
            raise ValueError("Strand does not have correct polarity.")

        self.template_polarity = [left_end, right_end]  # save just the polarity
        self.bases = template_strand[2:-2]  # save just the bases

        ## check the strand is DNA
        DNA, _ = self.DNA_or_RNA(self.bases)
        if not DNA:
            raise ValueError("Strand is not DNA")

    def DNA_or_RNA(self, bases):
        DNA_bases = set("ATGC")
        RNA_bases = set("AUGC")

        valid_DNA = set(bases) <= DNA_bases
        valid_RNA = set(bases) <= RNA_bases
        return valid_DNA, valid_RNA

    def DNA_to_RNA(self, bases):
        # T -> U
        return bases.replace("T", "U")

    def RNA_to_DNA(self, bases):
        # U -> T
        return bases.replace("U", "T")
