import json
import pysam
import pandas
from tqdm import tqdm
from pandas import ExcelWriter

known_sites = [
    241, 2790, 670, 2832, 3037, 4184, 4321, 5386, 8393, 9344, 9424,
    9534, 9866, 10029, 10198, 10447, 10449, 11537, 12880, 13195,
    14408, 15240, 15714, 17410, 18163, 19955, 20055, 21618, 21762,
    21846, 21987, 22200, 22578, 22673, 22674, 22679, 22686, 22688,
    22775, 22786, 22813, 22882, 22992, 22995, 23013, 23040, 23048,
    23055, 23063, 23075, 23202, 23403, 23525, 23599, 23604, 23854,
    23948, 24130, 24424, 24469, 24503, 25000, 25584, 26060, 26270,
    26530, 26577, 26709, 26858, 26259, 27382, 27383, 27384, 27807,
    28271, 28311, 28881, 28882, 28883, 29510,
]
known_files = pandas.read_csv(
    'file_name', delimiter='\t', encoding='utf-8', low_memory=False
)
writer = ExcelWriter('Output.xlsx')

for file in tqdm(known_files['files']):
    bamfile = pysam.AlignmentFile(file, "rb")
    df = pandas.DataFrame(index=['A', 'T', 'G', 'C', 'Total'], columns=known_sites)
    for site in known_sites:
        A, C, G, T = bamfile.count_coverage(
            bamfile.get_reference_name(0), start=site-1, stop=site, quality_threshold=0
        )
        df[site]['A'] = A[0]
        df[site]['T'] = T[0]
        df[site]['G'] = G[0]
        df[site]['C'] = C[0]
        df[site]['Total'] = A[0] + T[0] + G[0] + C[0]
    df.style.highlight_max().to_excel(writer, file.split('/')[1], header=True)

writer.save()
