import lzma
import warnings


def load_target_variants(loci_xz_path):
    loci = {}
    with lzma.open(loci_xz_path, 'rt') as infile:
        for line in infile:
            chrom, pos, ref, alt = line.rstrip().split('\t')
            loci[chrom + ':' + pos] = (ref, alt)
    return loci


def parse_vcf(vcf_stream, target_vars):
    for line in vcf_stream:
        if line[0] == '#':
            continue

        chrom, pos, _, ref, alts_string, _, _, _, fmt_fields, data_fields = line.rstrip().split('\t')

        locus = chrom + ':' + pos
        if locus not in target_vars:
            warnings.warn('Non-target locus {} encountered'.format(locus))
            continue
        locus_variants = target_vars[locus]

        if locus_variants[0] != ref:
            warnings.warn('Reference mismatch for locus {}'.format(locus))
            continue

        alts = alts_string.split(',')
        if locus_variants[1] not in alts:
            # Suppress output of loci without the target alt allele
            continue

        data = dict(zip(fmt_fields.split(':'), data_fields.split(':')))

        if 'AD' not in data:
            warnings.warn('Locus {} missing required AD FORMAT field'.format(locus))
            continue
        if 'PL' not in data:
            warnings.warn('Locus {} missing required PL FORMAT field'.format(locus))
            continue

        ad = data['AD'].split(',')
        ad_ref = int(ad[0])
        ad_alt = int(ad[alts.index(locus_variants[1])+1])

        # Keep only variants that are convincingly het.  We define this as
        # min(pl[0], pl[2:]) >= 20
        pl = [int(x) for x in data['PL'].split(',')]
        best_nonhet_pl = min(pl[0], min(pl[2:]))
        if best_nonhet_pl < 20:
            continue

        print('{}\t{}\t{}\t{}'.format(chrom, pos, ad_ref + ad_alt, ad_alt))


if __name__ == '__main__':
    import sys
    target_vars = load_target_variants(sys.argv[1])
    parse_vcf(sys.stdin, target_vars)

