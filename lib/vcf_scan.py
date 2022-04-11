import argparse
from utils import *
import re
import os

# Example usage: 
# python vcf_scan.py --wkdir data ## /home/nas210/missmi/Results/MMI001_00c4/variant_call \
#                    --id 001 ## 00c4 \
#                    --out vcf_scan_out \
#                    --report new_out ## /MMI001/MMI001_coverage_check_filter.txt
parser = argparse.ArgumentParser()
parser.add_argument('--wkdir', help='data directory', required=True)
parser.add_argument('--id',help='input subject id', required=True)
parser.add_argument('--out', help='output directory', required=True)
#parser.add_argument('--report', help='low coverage rport directory', required=True)
#parser.add_argument('--fast', help='run vcf scan without rescan', type=bool, default=False)
#parser.add_argument('--head', help="file head", required=True)

args = parser.parse_args()


subject_id = args.id
wk_dir = args.wkdir+'/variant_call'
out_dir = args.out+'/vcf_scan'
# v1
report_dir = args.out+'/'+subject_id+'_coverage_check.txt'
# v2
#report_dir = args.wkdir+'/MMI'+subject_id+'/analyze_acmg/MMI'+subject_id+'.low_coverage_candidate.txt'
#fast = args.fast

with open(report_dir,'r') as h:
    report_rows = h.readlines()

#if not fast:
if not os.path.isfile(f'{out_dir}/indel_only.recode.vcf'):
    cmd_list = [f"mkdir -p {out_dir}",
            f"vcftools --gzvcf {wk_dir}/{subject_id}.decomposed.normalized.vcf.gz --keep-only-indels --recode --out indel_only",
            f"mv indel_only.recode.vcf {out_dir}"]

    for cmd in cmd_list:
        call_cmd(cmd)

deletions = parse_vcf(f'{out_dir}/indel_only.recode.vcf', find_deletion)
print(f"Deletion number:  {len(deletions)}")
# build index for acceleration
vcf_idx = {}
for i,item in enumerate(deletions):
    k = item.split('\t')[0]
    if k not in vcf_idx.keys():
        vcf_idx[k] = i

# Search 'chr start end' from report
def search_region(line):
    '''
    check vcf dels if they overlapped the low coverage report 
    '''
    chr_, start, end = re.split(':|-|\s', line)[:3]
    key_list = list(vcf_idx.keys())
    chr_list_idx = key_list.index(chr_)
    chr_end = key_list[chr_list_idx+1]
    collect_vcf = []
    for item in deletions[vcf_idx[chr_]:vcf_idx[chr_end]]:
        del_info = item.split('\t')
        if del_info[0]==chr_:
            report_pos = (int(start),int(end))
            vcf_pos = (int(del_info[1]),int(del_info[1])+len(del_info[3])-1)
            if is_overlap(report_pos, vcf_pos):
                collect_vcf.append(item) 
    return collect_vcf

with open(out_dir+'/report_vcf_del.txt','w') as h:
    line_cnt, report_in = 0, 0
    for row in report_rows:
        results = search_region(row)
        if len(results)!=0:
            h.write(row)
            for res in results:
                h.write(f"{res}\n")
            line_cnt+=len(results)
            report_in+=1

# remove file if file is empty
#cmd = f"[ -s {out_dir}/report_vcf_del.txt ] || rm {out_dir}/report_vcf_del.txt"
#call_cmd(cmd)

print(f"Low coverage report number: {len(report_rows)}")
print(f"Low coverage report in vcf deletion: {report_in}")
print(f"vcf deletions in low coverage report: {line_cnt}")
print("vcf scan complete.")


