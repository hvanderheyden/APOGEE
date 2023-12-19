from glob import glob
import os
from collections import OrderedDict

# paf = /path/to/file.paf
paf_folder="/media/herve/10TB/Apogee/6_pipeline_mock/6_spagetthi/Scrubb_refilt/filteredPAFs"
out_file="/media/herve/10TB/Apogee/6_pipeline_mock/6_spagetthi/filtered_otu.tsv"


# Set minimum confidence
min_conf = 0


# List all paf files in folder
paf_list = glob(paf_folder + '/*.paf')

# Use dictinary to store date
# Key: target (identification)
# Value: list of mapping quality scores from reads matching that target
db_dict = dict()
sample_list = list()

for paf in paf_list:
    sample_name = os.path.basename(paf).split('-')[0]
    sample_list.append(sample_name)

    with open(paf, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            info_list = line.split('\t')

            # Extract columns of interest from paf file
            # query_name = info_list[0]
            target_name = info_list[5]
            map_qual = info_list[11]

            # Create dictionary entry if target seen for first time
            if target_name not in db_dict:
                db_dict[target_name] = dict()
                db_dict[target_name]['quals'] = [map_qual]  # Add Q-score to list
                db_dict[target_name]['n'] = 1  # Add to total mapping count
               
            else:
                db_dict[target_name]['quals'].append(map_qual)  # Add Q-score to list
                db_dict[target_name]['n'] += 1
               

            if sample_name not in db_dict[target_name]:
                db_dict[target_name][sample_name] = 1  # Add to sample count
            else:
                db_dict[target_name][sample_name] += 1

            # if sample_name not in sample_dict:
            #     sample_dict[sample_name] = dict()
            # if target_name not in sample_dict[sample_name]:
            #     sample_dict[sample_name][target_name] = 1  # Add to sample count
            # else:
            #     sample_dict[sample_name][target_name] += 1


def ave_prob(map_qual_list):
    # https://www.biostars.org/p/295932/
    if map_qual_list:
        return sum([10**(int(q) / -10) for q in map_qual_list]) / len(map_qual_list)
    else:
        return None


def prob_2_conf(prob):
    return round(100.0 - (prob * 100), 2)


# Create output file with confidense level
with open(out_file, 'w') as f:
    f.write('Accession\t{}\tTotalCount\tMappingConfidence\n'.format('\t'.join(sample_list)))  # File header
    for target_name, my_dict in db_dict.items():
        map_qual_list = my_dict['quals']
        n_reads = my_dict['n']

        # Compute confience
        prob = ave_prob(map_qual_list)
        conf = prob_2_conf(prob)

        # Get sample values for that specific target
        count_list = list()
        # for sample_name, count_dict in sample_dict.items():
        #     for t_name, count in count_dict.items():
        #         if t_name == target_name:
        #             count_list.append(count)
        for sample_name in sample_list:
            if sample_name in my_dict:
                count_list.append(my_dict[sample_name])
            else:
                count_list.append(0)

        # Write ID and conf value to file
        if int(conf) >= min_conf:
            f.write('{}\t{}\t{}\t{}\n'.format(target_name, '\t'.join([str(x) for x in count_list]), n_reads, conf))