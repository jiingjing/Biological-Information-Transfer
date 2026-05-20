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

    def base_pair_DNA(self, bases):
        # A - T
        # G - C
        base_table = str.maketrans("ATGC", "TACG")
        base_pair = bases.translate(base_table)
        return base_pair

    def base_pair_RNA(self, bases):
        # A - U
        # G - C
        base_table = str.maketrans("AUGC", "UACG")
        base_pair = bases.translate(base_table)
        return base_pair

    def match_codons(
        self, bases, RNA_polarity, name_type="Short Name", frame=1, terminus=True
    ):
        # split into triplets based on reading frame
        # 1 = default, start from first amino acid
        n = 3
        triplets = [(bases[i : i + n]) for i in range(frame - 1, len(bases) - frame, n)]

        protein = []
        for t in triplets:
            if RNA_polarity[1] == "5'":  # if 3'-5', read right to left i.e. backwards
                t = t[::-1]  # reverse the codon i.e. AUG -> GUA

            codon = self.codon_df.loc[self.codon_df["Codon"] == t, name_type].item()
            protein.append(codon)

        # check RNA_polarity and insert N/C terminus
        if terminus:
            if RNA_polarity[1] == "5'":  # if 3'-5', read right to left i.e. backwards
                protein.append("N")  # add N at end
                protein.insert(0, "C")  # add C at start
            else:  # if 5'-3', read left to right i.e. forwards
                protein.append("C")  # add C at end
                protein.insert(0, "N")  # add N at start

        return protein

    def replicate_DNA(self):
        replicated = self.base_pair_DNA(self.bases)
        return self.template_polarity[1] + replicated + self.template_polarity[0]

    def transcribe_RNA(self):
        transcribed = self.DNA_to_RNA(self.bases)
        transcribed = self.base_pair_RNA(transcribed)
        return self.template_polarity[1] + transcribed + self.template_polarity[0]

    def translate_RNA(
        self, rna_bases, RNA_polarity=None, name_type="Short Name", terminus=True
    ):
        if RNA_polarity is None:
            RNA_polarity = self.template_polarity[::-1]
        return self.match_codons(
            rna_bases, RNA_polarity, name_type=name_type, terminus=terminus
        )

    def info_transfer(self):
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

        # Translation
        proteins = self.translate_RNA(rna_bases=transcribed[2:-2])
        proteins = "-".join(proteins)
        print(f"Protein Translation: {proteins}")

        return replicated, transcribed, proteins
