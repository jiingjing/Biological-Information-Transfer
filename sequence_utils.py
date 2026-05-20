def gc_content(bases):
    """Return GC content as a percentage for a nucleotide sequence."""
    GC_count = bases.count("G") + bases.count("C")
    total_count = len(bases)
    if total_count == 0:
        return 0.0
    return (GC_count / total_count) * 100
