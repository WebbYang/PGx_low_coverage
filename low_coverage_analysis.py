import argparse
import os
from lib.utils import *
# call_cmd, low_coverage_check, annotate_clinvar_variant, count_ClinVar, annotate_exon, low_coverage_report
'''
args: bam_dir, subject_id, out_dir
example: python3 low_coverage_analysis.py --id GB0001 --out GB0001 --sex M

'''

parser = argparse.ArgumentParser()
parser.add_argument('--id',help='input subject id', required=True)
parser.add_argument('--out', help='output directory', required=True)
parser.add_argument('--sex', help='M of F', required=True)
parser.add_argument('--vd', help='vcf deletion  of reported region from called variant')

args = parser.parse_args()
subject_id = args.id
out_dir = args.out
sex = args.sex

# input file
#bam_dir = '/home/missmi/Missmi_Intelligence/Results/'
#bam_file = f'{bam_dir}/{subject_id}/variant_call/{subject_id}.recaled.bam' # realigned
bam_dir = '/home/dna/webb/low_coverage_test/open_bam'
bam_file = f'{bam_dir}/{subject_id}.realigned.bam'
missmi_lib_dir = '/home/missmi/bin/missmi.scripts/missmi_tools/missmi_tools/missmi_lib/workflows/missmi_parse_pgx/bin/lib'
pgx20_bed = f'{missmi_lib_dir}/db_pgx20/20210501/GRCh37_latest_genomic.gene.pgx20_latest.rename.slopped.bed'

# output files
pgx_bam_file = f'{out_dir}/{subject_id}_pgx.bam'
low_coverage_candidate = f'{out_dir}/{subject_id}.low_coverage_candidate.txt'
annotate_clinvar_file = f'{out_dir}/{subject_id}_clinvar'
annotate_exon_file = f'{out_dir}/{subject_id}_exon'
low_coverage_file = f'{out_dir}/{subject_id}_coverage_check.txt'

# Extract bam
extract = "/home/missmi/bin/extract_pgx_bam"
def extract_pgx():
    print("Start extracting PGx region.")  
    cmd_list = [f"mkdir -p {out_dir}",
                f'{extract} {pgx_bam_file} {bam_file}']
    for cmd in cmd_list:
        _ = call_cmd(cmd)
    
    print(f'Extract {bam_file} to {out_dir}/{subject_id}_pgx.bam done!')

# Get low-coverage region
def get_low_coverage():
    res = low_coverage_check(pgx_bam_file, pgx20_bed)
    with open(low_coverage_candidate, 'w') as wh:
        for line in res:
            start, end = line.split(':')[1].split('-')
            dist = int(end) - int(start) + 1
            wh.write(f'{line}\t{dist}\n')

    print(f'Generate {low_coverage_candidate} done!')

def annotate_region():
    print('Start annotating ClinVar variants')   
    _ = annotate_clinvar_variant(low_coverage_candidate, annotate_clinvar_file)
    print('Start annotating exon region') 
    _ = annotate_exon(low_coverage_candidate, annotate_exon_file)

def report(filt=True):
    print('Start writing the report.')
    low_coverage_report(low_coverage_candidate, annotate_clinvar_file, annotate_exon_file, low_coverage_file)
    if filt:
        filtered_file = report_filter()
        return count_ClinVar(annotate_clinvar_file, filtered_file)
    return count_ClinVar(annotate_clinvar_file, low_coverage_file)

def report_filter():
    likely_deletion_check(pgx_bam_file, out_dir, subject_id)
    chrX_check(sex, pgx_bam_file, out_dir, subject_id)
    if args.vd:
        del_file = args.vd
        drop_del(out_dir, subject_id, del_file)
    else:
        del_file = f"/home/dna/webb/low_coverage_test/low_coverage_out/{subject_id}/vcf_scan/report_vcf_del.txt"
        if os.path.isfile(del_file):
            drop_del(out_dir, subject_id, del_file)
        else:
            drop_del(out_dir, subject_id)
    return low_coverage_file.replace('.txt', '_filter.txt')

def main():
    extract_pgx()
    get_low_coverage()
    annotate_region()
    report()

if __name__ == "__main__":
    main()
