# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 01:15:35 2023

@author: Bhumika
"""

import os
import re
import glob
import pandas as pd

#----------satya---------------#
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to process data with specific directory path')
    parser.add_argument('directory', type=str, help='Working directory path')
    args = parser.parse_args()

    os.chdir(args.directory)

    print(f"Changed working directory to: {os.getcwd()}")


#sr fge file
file_path = next(glob.iglob("*_SR_fge.txt"), None)

if file_path:
    SR_fge = pd.read_csv(file_path, delimiter="\t")
else:
    print('NO SR FGE file found')

#----------satya---------------#



#filtering genes having read-support >= 3        
rs_3 = SR_fge[SR_fge['supportRead'] >= 3].copy()

#filtering genes based on 11 genes ('ALK', 'ROS1', 'RET', 'MET','NTRK1','NTRK2','NTRK3','NRG1','NRG2','FGFR1','FGFR2','FGFR3')
fus_genes = ['ALK', 'ROS1', 'RET', 'MET','NTRK1','NTRK2','NTRK3','NRG1','NRG2','FGFR1','FGFR2','FGFR3']
panel_fus_gene_df1 = rs_3[rs_3['fusionName'].str.contains('|'.join(fus_genes))]

#reading the bedpe file in the current working directory
bedpe = pd.read_csv('splitReadInfo.bedpe', sep='\t')

#reading the backend db csv file in the current working directory
backend_db = pd.read_csv("https://raw.githubusercontent.com/atomikkus/FWES-DEM/main/Database/Exon_Intron_all_Region_FuSEQ.csv")

#--------------------------starting on bedpe file-------------------------------

# Splitting the name column in bedpe file by "--" and put the two parts in separate columns
bedpe[['name_split1', 'name_split2']] = bedpe['name'].str.split('--', n=1, expand=True)

# Loop over the specified fus_genes
filtered_data = pd.DataFrame()  # initialize an empty DataFrame
for fus_gene in fus_genes:
    # filtering the rows based on the name_split1 and name_split2 column if the specified fus_genes present and based on read_support >= 3
    filtered_data_fus_genes = bedpe[((bedpe['name_split1'] == fus_gene) | (bedpe['name_split2'] == fus_gene)) & (bedpe['score'] >= 3)]

    # appending the filtered data to the new DataFrame
    filtered_data = pd.concat([filtered_data, filtered_data_fus_genes])

# Dropping duplicates based on the 'name', 'chr1', 'start1', and 'end1' columns in bedpe file
unique_data_bedpe = filtered_data.drop_duplicates(subset=['name', 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'score'])

# Create an empty DataFrame to store only the required results from bedpe
bedpe_updated = pd.DataFrame(columns=['chrom1', 'end1', 'gene1', 'chrom2', 'start2', 'gene2', 'Fusion_name', 'Score/Read support'])

# Iterating over the unique_data_bedpe rows
for index, row in unique_data_bedpe.iterrows():
    # Split the name into two parts
    name_split = row['name'].split('--')
    name_split1 = name_split[0]
    name_split2 = name_split[1] if len(name_split) > 1 else None
    
    # Checking if name_split1 has any of the specified fus_genes
    if any(x in name_split1 for x in fus_genes):
        # Add the corresponding chrom, start, and end to the new_df 
        new_row = {'chrom1': row['chrom1'], 'start1': row['start1'], 'end1': row['end1'], 'gene1': name_split1, 'chrom2': row['chrom2'], 'start2': row['start2'], 'end2': row['end2'], 'gene2': name_split2, 'Fusion': 'NA', 'Fusion_name': row['name'], 'Score/Read support': row['score']}
        new_row_df = pd.DataFrame([new_row])  # Create a DataFrame from the dictionary
        bedpe_updated = pd.concat([bedpe_updated, new_row_df], ignore_index=True)
    
    # Check if name_split2 has any of the specified fus_genes
    if name_split2 and any(x in name_split2 for x in fus_genes):
        # Add the corresponding chrom, start, and end to the new_df
        new_row = {'chrom1': row['chrom2'], 'start1': row['start2'], 'end1': row['end2'], 'gene1': name_split2, 'chrom2': row['chrom1'], 'start2': row['start1'], 'end2': row['end1'], 'gene2': name_split1, 'Fusion': row['name'], 'Fusion_name': f"{row['name_split2']}--{row['name_split1']}", 'Score/Read support': row['score']}
        new_row_df = pd.DataFrame([new_row])  # Create a DataFrame from the dictionary
        bedpe_updated = pd.concat([bedpe_updated, new_row_df], ignore_index=True)


#-------------------------------mapping filtered_bed.pe file with fge file-------------------------------------
map_bedpe_updated_fge = pd.DataFrame()
for index, row in bedpe_updated.iterrows():
    match_count = 0  # Counter variable to track the number of matches
    for index1, row1 in panel_fus_gene_df1.iterrows():
        if row["Fusion_name"] == row1["fusionName"]:
            new_column = {"chrom1": row["chrom1"], "start1": row["start1"], "end1": row["end1"], "gene1": row["gene1"], "chrom2": row["chrom2"], "start2": row["start2"], "end2": row["end2"], "gene2": row["gene2"], "Fusion_name": row["Fusion_name"], "Fusion": row["Fusion"], "supportRead": row['Score/Read support'], "Read_name": row1['header']}
            new_column_df = pd.DataFrame([new_column])  # Create a DataFrame from the dictionary
            map_bedpe_updated_fge = pd.concat([map_bedpe_updated_fge, new_column_df], ignore_index=True)
            match_count += 1
            
        elif row["Fusion"] == row1["fusionName"]:
            new_column = {"chrom1": row["chrom1"], "start1": row["start1"], "end1": row["end1"], "gene1": row["gene1"], "chrom2": row["chrom2"], "start2": row["start2"], "end2": row["end2"], "gene2": row["gene2"], "Fusion_name": row["Fusion_name"], "Fusion": row["Fusion"], "supportRead": row['Score/Read support'], "Read_name": row1['header']}
            new_column_df = pd.DataFrame([new_column])  # Create a DataFrame from the dictionary
            map_bedpe_updated_fge = pd.concat([map_bedpe_updated_fge, new_column_df], ignore_index=True)
            match_count += 1

    # If match count is more than 1, a new row will be created
    if match_count > 1:
        new_column = {"chrom1": row["chrom1"], "start1": row["start1"], "end1": row["end1"], "gene1": row["gene1"], "chrom2": row["chrom2"], "start2": row["start2"], "end2": row["end2"], "gene2": row["gene2"], "Fusion_name": row["Fusion_name"], "Fusion": row["Fusion"], "supportRead": row['Score/Read support'], "Read_name": row1['header']}
        new_column_df = pd.DataFrame([new_column])  # Create a DataFrame from the dictionary
        map_bedpe_updated_fge = pd.concat([map_bedpe_updated_fge, new_column_df], ignore_index=True)


# Create an empty DataFrame to store the final results
map_with_backend_db = pd.DataFrame(columns=['#FusionGene', 'ReadSupport', 'LeftBreakpoint', 'gene1', 'RightBreakpoint', 'gene2', 'ReadNames', 'Fusion', 'Exon_Num', 'Exon_Cat'])


#--------------------mapping of map_bedpe_updated_fge file to backend_db csv file-------------------------------------
for index, row in map_bedpe_updated_fge.iterrows():
    for index1, row1 in backend_db.iterrows():
        if row["gene1"] == row1["Gene"]:
            # Check if the chromosome is the same in both the files
            if row['chrom1'] == row1['Chrom']:
                # Check if the start or end location in bedpe file is within the range of locations in backend database
                if (row1['Exon_Start'] >= row['end1'] >= row1['Exon_End']) or (row1['Exon_Start'] <= row['end1'] <= row1['Exon_End']):
                    if row['Fusion'] == "NA":
                        new_row = {'#FusionGene': row['Fusion_name'], 'ReadSupport': row['supportRead'], 'LeftBreakpoint': f"{row['chrom1']}:{row['end1']}", 'gene1': row['gene1'], 'RightBreakpoint': f"{row['chrom2']}:{row['start2']}", 'gene2': row['gene2'] ,'ReadNames': row['Read_name'], 'Fusion': row['Fusion'], 'Exon_Num': row1['Exon_Num'],'Exon_Cat': row1['Exon_Cat']}
                        new_row_df = pd.DataFrame([new_row])  # Create a DataFrame from the dictionary
                        map_with_backend_db = pd.concat([map_with_backend_db, new_row_df], ignore_index=True)
                    else:
                        new_row = {'#FusionGene': row['Fusion'],'ReadSupport': row['supportRead'], 'LeftBreakpoint': f"{row['chrom2']}:{row['end2']}", 'gene1': row['gene2'], 'RightBreakpoint': f"{row['chrom1']}:{row['start1']}", 'gene2': row['gene1'] ,'ReadNames': row['Read_name'], 'Fusion': row['Fusion'], 'Exon_Num': row1['Exon_Num'],'Exon_Cat': row1['Exon_Cat']}
                        new_row_df = pd.DataFrame([new_row])  # Create a DataFrame from the dictionary
                        map_with_backend_db = pd.concat([map_with_backend_db, new_row_df], ignore_index=True)
                else:
                    pass


# To group by the specified columns and concatenate values in Exon_Cat and Exon_Num columns
df_merged = map_with_backend_db.groupby(['#FusionGene', 'ReadSupport', 'LeftBreakpoint', 'gene1', 'RightBreakpoint', 'gene2', 'ReadNames', 'Fusion']).agg({'Exon_Cat': lambda x: ', '.join(x.astype(str)),'Exon_Num': lambda x: ', '.join(x.astype(str))}).reset_index()


#---------satya---------------------------Filter based on exon numbers--------------------------------------

# Define a dictionary for gene-specific exon numbers
exon_numbers = {
    'ALK': ['17', '18', '19', '20', '21', '22'],
    'RET': ['2','3','4','5','6','7','8','9','10','11','12','13','14'],
    'MET': ['12','13','14','15','16'],
    'ROS1': ['29','30','31','32','33','34','35','36','37'],
    'NTRK1': ['2','3','4','5','6','7','8','9','10','11','12','13','14'],
    'NTRK2': ['8','9','10','11','14','15','16','17'],
    'NTRK3': ['3','4','5','12','13','14','15','16'],
    'FGFR1': ['6','7','8','13','14','15','16','17','18'],
    'FGFR2': ['14','15','16','17','18','19'],
    'FGFR3': ['13','14','15','16','17','18','19'],
    'NRG1': ['1','2','3','4','5','6','7'],
    'NRG2': ['1','2','3','4','5','6','7']
}

# Filter based on gene-specific exon numbers
filtered_exons = pd.DataFrame(columns=df_merged.columns)

for index, row in df_merged.iterrows():
    fusiongene = row['#FusionGene']
    exon_cat = row['Exon_Cat']
    
    for gene, exons in exon_numbers.items():
        if gene in fusiongene and any(x in exon_cat for x in exons):
            row_df = pd.DataFrame([row])
            filtered_exons = pd.concat([filtered_exons, row_df], ignore_index=True)
            break  # Break once a matching gene and exon combination is found

#reading the known_fusions csv file in the current working directory
known_Fusions = pd.read_csv("https://raw.githubusercontent.com/atomikkus/FWES-DEM/main/Database/known_Fusions.csv")

#storing formats to be sent to reporting team
final_df = pd.DataFrame(columns=['#FusionGene', 'ReadSupport', 'LeftBreakpoint', 'RightBreakpoint', 'ReadNames'])

#format to check the FusionStatus ("Unknown/known")
final_df_exons = pd.DataFrame(columns=['#FusionGene', 'ReadSupport', 'LeftBreakpoint', 'RightBreakpoint', 'ReadNames', 'Exons', 'FusionStatus'])


for index, row in filtered_exons.iterrows():
    fusion_matched = False  # Flag to check if a match is found
    
    for index1, row1 in known_Fusions.iterrows():
        if row["#FusionGene"] == row1["KnownFusions"]:
            new_row2 = {'#FusionGene': row['#FusionGene'],'ReadSupport': row['ReadSupport'],'LeftBreakpoint': row['LeftBreakpoint'],'RightBreakpoint': row['RightBreakpoint'],'ReadNames': row['ReadNames'], 'Exons':row['Exon_Cat'], 'FusionStatus': 'KnownFusions'}
            new_row2_df2 = pd.DataFrame([new_row2])  
            final_df_exons = pd.concat([final_df_exons, new_row2_df2], ignore_index=True)
            
            fusion_matched = True  # Setting the flag to True to indicate a match
            break  # Break out of the inner loop

    # Adding rows to final_df and final_df_exons for the case where no match was found
    new_row1 = {'#FusionGene': row['#FusionGene'],'ReadSupport': row['ReadSupport'],'LeftBreakpoint': row['LeftBreakpoint'],'RightBreakpoint': row['RightBreakpoint'],'ReadNames': row['ReadNames']}
    new_row1_df = pd.DataFrame([new_row1])  
    final_df = pd.concat([final_df, new_row1_df], ignore_index=True)
    
    if not fusion_matched:
        new_row2 = {'#FusionGene': row['#FusionGene'],'ReadSupport': row['ReadSupport'],'LeftBreakpoint': row['LeftBreakpoint'],'RightBreakpoint': row['RightBreakpoint'],'ReadNames': row['ReadNames'], 'Exons':row['Exon_Cat'], 'FusionStatus': 'Unknown'}
        new_row2_df2 = pd.DataFrame([new_row2])  
        final_df_exons = pd.concat([final_df_exons, new_row2_df2], ignore_index=True)
            
#--------------------------------------------------------------------------------

# To remove the trailing ".0" only if there are more than one trailing zeros from the output
def remove_trailing_zeros(value):
    if value.endswith('.0'):
        value = value.rstrip('0').rstrip('.')
    return value

# To remove the trailing ".0" from 'LeftBreakpoint' and 'RightBreakpoint' columns
final_df['LeftBreakpoint'] = final_df['LeftBreakpoint'].apply(remove_trailing_zeros)
final_df['RightBreakpoint'] = final_df['RightBreakpoint'].apply(remove_trailing_zeros)

final_df_exons['LeftBreakpoint'] = final_df_exons['LeftBreakpoint'].apply(remove_trailing_zeros)
final_df_exons['RightBreakpoint'] = final_df_exons['RightBreakpoint'].apply(remove_trailing_zeros)

#--------------------------------------------------------------------------------

#-----------------------------------------Generating sample wise final output files--------------------------------
#getting current working directory
path = os.getcwd()

# Extracting the sample_id from the working directory
file_name = os.path.splitext(os.path.basename(path))[0]

# Extracting the desired string (Sample_ID as named in the current working directory)
desired_string = re.search(r'(.+?)_fusions', file_name).group(1)

##########################  adding sample_id to the fuseq-output.xlsx file  #################################
#inserting Sample_ID column to the final output
map_with_backend_db.insert(0, 'Sample_ID', '')

#inserting the sample_ID to the respective fusions
map_with_backend_db['Sample_ID'] = desired_string

# Output of all fusions having read support >= 3
# map_with_backend_db.to_excel(f'{desired_string}_fuseq-output.xlsx', index=False)

##########################  sample adding done  #################################

#output for bioinfo validation
#filtered_rs.to_excel(f'{desired_string}_DNA_fuseq-output.xlsx', index=False)

#output for reporting team
final_df.to_excel(f'{desired_string}_DNA_FUS_P.xlsx', index=False)

#output for FusionStatus ("Unknown/known")
final_df_exons.to_excel(f'{desired_string}_DNA_FUS_Exons_P.xlsx', index=False)

#-------------------------------------------------satya-------------------------------------------#


# Reading the exon coordinates file into a DataFrame
exon_df = backend_db.copy()

# Reading the output file into a DataFrame (assuming the output is stored in a file named 'output.csv')
output_df = final_df_exons.copy()

# Function to adjust coordinates based on 'Exons' column using 'exon_coord.csv'
def adjust_coordinates(row):
    if 'Intron' in row['Exons']:
        exon_numbers = [int(s.split('.')[0]) for s in row['Exons'].split('_')[-1].split('-')]
        
        # Fetching exon coordinates from exon_df for the relevant gene
        gene_exons = exon_df[(exon_df['Gene'] == row['Gene']) & (exon_df['Exon_Num'].isin(exon_numbers))]
        
        # Checking the distance to Exon_Start and Exon_End for each relevant exon
        gene_exons['Dist_Start'] = abs(gene_exons['Exon_Start'] - row['Gene_Start'])
        gene_exons['Dist_End'] = abs(gene_exons['Exon_End'] - row['Gene_End'])
        
        # Choosing the exon with the nearest start or end depending on the minimum distance
        nearest_exon = gene_exons.loc[
            (gene_exons['Dist_Start'] == gene_exons['Dist_Start'].min()) |
            (gene_exons['Dist_End'] == gene_exons['Dist_End'].min())
        ]
        
        return nearest_exon.iloc[0]['Exon_Start'] if nearest_exon.iloc[0]['Dist_Start'] < nearest_exon.iloc[0]['Dist_End'] else nearest_exon.iloc[0]['Exon_End']
    else:
        return row['Gene_Start']

# Apply the function to create a new column 'Adjusted_Coordinates'
output_df[['Gene1','Gene2']] = output_df['#FusionGene'].str.split('--', expand=True)
output_df[['chrL', 'GC1']] = output_df['LeftBreakpoint'].str.split(':', expand=True)
output_df['GC1'] = output_df['GC1'].astype(int)


output_df['Adjusted_Coordinates'] = output_df.apply(adjust_coordinates, axis=1)
output_df['Adjusted_Coordinates'] = output_df['chrL'] + ':' + output_df['GC1'].astype('str')


# Save the updated DataFrame to a new file
output_df.drop(['chrL','GC1','Gene1','Gene2'], axis=1, inplace=True)
output_df.to_csv('output_with_adjusted_coordinates.csv', index=False)
