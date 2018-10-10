import sys
# This program is intended as an example of solving the final project in
# Introduction to Programming part of SIRE504.
# Written by Harald Grove, 2017
# e-mail: haraldgrove@gmail.com
# The program is intended to only use basic python functions and only uses the
# sys-module for access to the command line via the sys.argv list variable.

def explain_use():
    """
    Prints a short explanation to the user about how to run the program.
    """
    print('Usage: parseVCF command vcf-file [options]')
    print(' Commands: stats tag translate')

def explain_command(command):
    """
    Prints a short task specific help text, including expected input parameters
    """
    if command == 'stats':
        print('Usage: parseVCF', command, 'vcf-file')
    elif command == 'tag':
        print('Usage: parseVCF', command, 'vcf-file tag')
    elif command == 'translate':
        print('Usage: parseVCF', command, 'vcf-file output-prefix')
    else:
        print('Unknown command', command)
    
def gather_statistics(vcffile):
    """ 
    Invoked with the 'stats' option.
    Prints the number of samples in the file, together with the number of SNPS and INDELs.
    Will also report if any TAGs used are not described in the header section.
    """
    with open(vcffile, 'r') as fin:
        db = {}
        bad_terms = []
        sample_names = []
        indels = 0
        snps = 0
        for line in fin:
            if line.startswith('##INFO') or line.startswith('##FORMAT'):
                ind = line.index('ID=')
                key = line[ind+3:ind+5]
                db[key] = 1
                continue
            if line.startswith('#CHROM'):
                l = line.strip().split()
                sample_names = line.strip().split()[9:]
                if l[1] != 'POS' or l[2] != 'ID' or l[3] != 'REF' or l[4] != 'ALT' or\
                   l[5] != 'QUAL' or l[6] != 'FILTER' or l[7] != 'INFO':
                    print('Unexpected name in the column header')
                    print('{}'.format('\t'.join(l[0:10])))
            if line.startswith('#'):
                continue
            chrom, pos, id_, ref, alt, qual, filter, info, format_, *samples = line.strip().split()
            info_l = info.split(';')
            format_l = format_.split(':')
            for element in info_l+format_l:
                if element[:2] not in db and element[:2] not in bad_terms:
                    bad_terms.append(element[:2])
                if element == 'VT=SNP':
                    snps += 1
                elif element == 'VT=INDEL':
                    indels += 1
    print('Samples:', len(sample_names))
    print('SNPs:', snps)
    print('INDELs:', indels)
    if len(bad_terms) > 0:
        print('Tags missing from header section:', bad_terms)

def find_term(vcffile, term):
    """
    Extract information from one TAG for each variant and displays this, together
    with the chromosome number and position.
    The code also handles the case when the TAG is a FORMAT TAG, even if this was not asked for.
    """
    with open(vcffile, 'r') as fin:
        tag_type = ''
        notfound = False
        print('CHROM\tPOS',term)
        for line in fin:
            if line.startswith('#'):
                if '<ID='+term in line and line.startswith('##INFO'):
                    tag_type = 'INFO'
                elif '<ID='+term in line and line.startswith('##FORMAT'):
                    tag_type = 'FORMAT'
                continue
            if tag_type == '':
                print('Tag',term,'was not found in either INFO or FORMAT.')
                return
            chrom, pos, id_, ref, alt, qual, filter, info, format_, *samples = line.strip().split()
            if tag_type == 'INFO':
                info_l = info.split(';')
                for il in info_l:
                    if term in il:
                        print(chrom, pos, il[3:])
                        break
                else:
                    print(chrom, pos, '.')
            elif tag_type == 'FORMAT':
                format_l = format_.split(':')
                ind = format_l.index(term)
                output = [chrom, pos]
                for sample in samples:
                    sl = sample.split(':')[ind]
                    output.append(sl)
                print('\t'.join(output))
            else:
                notfound = True
        if notfound:
            print('Tag',term,' was missing for one or more variants.')
    
def translate_plink(vcffile, prefix):
    """
    Extracts genotypes for each SNP for all samples and converts it to a simplified
    PLINK output.
    The output also includes a second file displaying the positional information
    for each SNP.
    """
    mapfile = prefix+'.map'
    pedfile = prefix+'.ped'
    sample_names = []
    with open(vcffile, 'r') as fin, open(mapfile, 'w') as fout:
        db = {}
        for line in fin:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                sample_names = line.strip().split()[9:]
                for name in sample_names:
                    db[name] = [name]
                continue
            chrom, pos, id_, ref, alt, qual, filter, info, format_, *samples = line.strip().split()
            if 'VT=SNP' not in info or len(ref+alt) != 2:
                continue
            fout.write('{}\t{}\t{}\n'.format(chrom, id_, pos))
            trans = {'0':ref, '1':alt, '.':'0'}
            for i, name in enumerate(sample_names):
                gt1,gt2 = samples[i][0],samples[i][2]
                db[name].append(trans[gt1])
                db[name].append(trans[gt2])
    with open(pedfile, 'w') as fout:
        for name, geno in db.items():
            fout.write('{}\n'.format('\t'.join(geno)))
            

def main():
    """
    Interpretes the command line, running the appropriate section depending on input.
    Also provides some helpful instructions if input is missing or does not match
    expectations.
    """
    if len(sys.argv) == 1:
        explain_use()
    elif len(sys.argv) == 2:
        explain_command(sys.argv[1])
    elif sys.argv[1] == 'stats':
        gather_statistics(sys.argv[2])
    elif sys.argv[1] == 'tag':
        if len(sys.argv) == 4:
            find_term(sys.argv[2], sys.argv[3])
        else:
            print('Incorrect number of options.')
    elif sys.argv[1] == 'translate':
        if len(sys.argv) == 4:
            translate_plink(sys.argv[2], sys.argv[3])
        else:
            print('Incorrect number of options.')
    else:
        print('Wrong input, please check your file-names and options and try again.')
        explain_use()

if __name__ == '__main__':
    main()
