import gzip
import json
import sys

def read_json_to_dict(dictfile: str) -> dict:
    with open(dictfile) as f:
        pos_count_dict = json.load(f)

    return pos_count_dict

def create_dict(filename: str, pos_count_dict: dict) -> dict:

    positions_with_patients = {}

    with gzip.open(filename, mode="rt") as file:

        for i, line in enumerate(file):
            #if i > 1000000:
            #    break

            sys.stdout.write("Progress: working on line %d   \r" % (i) )
            sys.stdout.flush()

            splitline = line.split("\t")

            if i == 0:
                header_list = splitline

                try:
                    normal_gt_index = header_list.index("normal_genotype")
                    tumor_gt_index = header_list.index("tumor_genotype")
                    pid_index = header_list.index("PID")
                    alt_index = header_list.index("ALT")
                    sample_control_index = header_list.index("sample_control")
                    sample_tumor_index = header_list.index("sample_tumor")
                    start_index = header_list.index("start")
                    end_index = header_list.index("end")
                    quality_score_index = header_list.index("quality_score")
                    reads_normal_index = header_list.index("reads_normal")
                    reads_tumor_index = header_list.index("reads_tumor")

                except ValueError:
                    print("element not in header")

            normal_gt = splitline[normal_gt_index]
            tumor_gt = splitline[tumor_gt_index]
            pid = splitline[pid_index]
            chrom = splitline[0]
            pos = splitline[1]
            key = chrom + "-" + pos

            if normal_gt == "1,1" and "0" in tumor_gt and key in pos_count_dict.keys():
                positions_with_patients[key] = {"count": pos_count_dict[key]}
                positions_with_patients[key][pid] = {
                    "ALT": splitline[alt_index],
                    "sample_control": splitline[sample_control_index],
                    "sample_tumor": splitline[sample_tumor_index],
                    "start": splitline[start_index],
                    "end": splitline[end_index],
                    "normal_genotype": normal_gt,
                    "tumor_genotype": tumor_gt,
                    "quality_score": splitline[quality_score_index],
                    "reads_normal": splitline[reads_normal_index],
                    "reads_tumor": splitline[reads_tumor_index]
                }

    return positions_with_patients

def main(filename: str, outfile: str, dictfile: str):

    pos_count_dict = read_json_to_dict(dictfile)
    positions_with_patients = create_dict(filename, pos_count_dict)
    with open(outfile, "w") as file:
        json.dump(positions_with_patients, file)

if __name__ == "__main__":

    filename = "/home/junkhann/daten/mmml_onek1k_all_patients.vcf.gz"
    outfile = "/home/junkhann/daten/positions_with_patients.json"
    dictfile = "/home/junkhann/daten/position_count.json"

    main(filename, outfile, dictfile)