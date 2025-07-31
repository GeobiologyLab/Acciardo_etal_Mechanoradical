import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib as mpl



ros = pd.read_csv('Data/ros_ko_heatmap_table_reordered.csv',sep=',',header=0)



myOrder = ['MAG_1','MAG_2','MAG_3','MAG_4','MAG_5','MAG_6','MAG_7','MAG_8','MAG_9','MAG_10', 'MAG_11','MAG_12','MAG_13','MAG_14']

taxonomy = [
    'Desulforudis','Hydrogenophaga','Hydrogenophaga', 'Alphaproteobacteria', 'Thermodesulfobacteriota',
   'Bacillota','Spirochaetota','Thermodesulfovibrionia',
    'Thermodesulfovibrionia','Nitrospirota','Actinomycetota','Thermodesulfobacteriota',
    'Alphaproteobacteria','Spirochaetota'
]

colors = {'Desulforudis':'a5cee0','Hydrogenophaga':'6a3e98','Thermodesulfobacteriota':'f58120','Alphaproteobacteria':'f69799','Thermodesulfovibrionia':'b7d885','Bacillota':'1f78b4','Spirochaetota':'808080','Nitrospirota':'808080','Actinomycetota':'808080'}


dataframes = []

for i, mag in enumerate(myOrder):
    id = mag
    input = 'MAGS/%s.kegg.list' % mag
    ros_row = ros[ros['genome'] == id]
    ros_row = ros_row.drop(columns=['genome'])
    print(input)
    if i == 0:
        df = pd.read_csv(input, sep="\t", header=0)
        ros_row.index = df.index
        merged_df = pd.concat([df, ros_row], axis=1)
    else:
        df = pd.read_csv(input, sep="\t", header=0)
        ros_row.index = df.index
        df = pd.concat([df, ros_row], axis=1)
        merged_df = pd.concat([merged_df, df], axis=0, ignore_index=True)
     
print(merged_df)

merged_df = merged_df.apply(pd.to_numeric, errors='coerce', axis=0)

merged_df['Glycoside Hydrolase*'] = ((merged_df['D-galacturonate epimerase'] ==1) | (merged_df['beta-N-acetylhexosaminidase'] ==1)).astype(int)
print(merged_df.columns.to_list())

merged_df['C-P Lyase*'] = ((merged_df['C-P lyase cleavage PhnJ'] ==1) | (merged_df['CP-lyase complex'] ==1)).astype(int)




mySubset = ['Cytochrome c oxidase', 'Ubiquinol-cytochrome c reductase','Cytochrome c oxidase, cbb3-type', 'Cytochrome bd complex',  'NAD-reducing hydrogenase', 
            'NiFe hydrogenase Hyd-1', 'CBB Cycle', 'Wood-Ljungdahl', 'Glycoside Hydrolase*',  'dissimilatory sulfate < > APS', 'dissimilatory sulfite < > APS',  'dissimilatory sulfite < > sulfide', 'thiosulfate oxidation', 'thiosulfate/polysulfide reductase',
            'sulfur dioxygenase',  'sulfide oxidation','nitrogen fixation', 'dissim nitrate reduction', 'DNRA', 'hydroxylamine oxidation', 'nitrite oxidation', 'polyhydroxybutyrate synthesis',
             'bidirectional polyphosphate', 'C-P Lyase*','Flagellum', 'Chemotaxis','Biofilm PGA Synthesis protein']
mySubset = ros.columns[1:].to_list() + mySubset


blues = ros.columns[1:].to_list() + ['Cytochrome c oxidase', 'Ubiquinol-cytochrome c reductase', 'Cytochrome c oxidase, cbb3-type', 'Cytochrome bd complex', 'NAD-reducing hydrogenase', 'NiFe hydrogenase Hyd-1']
greys =  ['Wood-Ljungdahl', 'Glycoside Hydrolase*']
oranges = ['dissimilatory sulfate < > APS', 'dissimilatory sulfite < > APS', 'dissimilatory sulfite < > sulfide', 'thiosulfate oxidation', 'thiosulfate/polysulfide reductase', 'sulfur dioxygenase', 'sulfide oxidation']
purples = ['nitrogen fixation', 'dissim nitrate reduction', 'DNRA', 'hydroxylamine oxidation', 'nitrite oxidation']
greens = ['polyhydroxybutyrate synthesis', 'bidirectional polyphosphate', 'C-P Lyase*', 'Flagellum', 'Chemotaxis', 'Biofilm PGA Synthesis protein']

all_groups = blues + greys + oranges + purples + greens
data = merged_df[all_groups]

fig, ax = plt.subplots(figsize=(len(all_groups) * 0.5, 8))

mask = np.full(data.shape, True)

groups = {
    'Blues': (blues, 'O & H'),
    'Greys': (greys, 'C'),
    'Reds': (oranges, 'S'),
    'Purples': (purples, 'N'),
    'Greens': (greens, 'Other')
}

norm = mpl.colors.Normalize(vmin=np.nanmin(data.values), vmax=np.nanmax(data.values))
cbar_axes = []

colorbar_positions = [0.82, 0.85, 0.88, 0.91, 0.94]  # more tightly packed
cbar_idx = 0

for cmap_name, (columns, legend_title) in groups.items():
    group_mask = np.full(data.shape, True)
    cols_idx = [data.columns.get_loc(c) for c in columns]
    group_mask[:, cols_idx] = False

    cmap = sns.color_palette(cmap_name, as_cmap=True)

    sns.heatmap(
        data,
        mask=group_mask,
        cmap=cmap,
        cbar=False,
        linewidths=0.5,
        linecolor='gray',
        square=True,
        ax=ax
    )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar_ax = fig.add_axes([colorbar_positions[cbar_idx], 0.3, 0.015, 0.4])  


    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.ax.set_title(legend_title, rotation=0, fontsize=10, pad=10)
    if legend_title != 'Other':
        cbar.ax.set_yticklabels([])  
        cbar.ax.tick_params(size=0)  



title_case_labels = [label.title() for label in all_groups]
for i,label in enumerate(title_case_labels):
    if label =='Dnra':
        title_case_labels[i] = 'DNRA'
    if ' C ' in label:
        label = label.replace(' C ', ' c ')
        title_case_labels[i] = label
    if 'Cbb' in label:
        title_case_labels[i] = 'CBB Cycle'
    if 'Pga' in label:
        label = label.replace('Pga', 'PGA')
        title_case_labels[i] = label
    if 'Bc' in label:
        label = label.replace('Bc', 'bc')
        title_case_labels[i] = label
    if 'And/Or' in label:
        label = label.replace('And/Or', 'and/or')
        title_case_labels[i] = label

ax.set_xticks(np.arange(len(all_groups)) + 0.5)
ax.set_xticklabels(title_case_labels, rotation=90, fontsize=12)


yticklabels = ax.get_yticklabels()

new_labels = []
for label in yticklabels:
    try:
        new_value = int(label.get_text()) + 1
        new_labels.append(str(new_value))
    except ValueError:
        new_labels.append(label.get_text())

ax.set_yticklabels(new_labels, fontsize=12)

for label, taxon in zip(ax.get_yticklabels(), taxonomy):
    if taxon in colors:
        label.set_color('#' + colors[taxon])  # Add "#" to hex code
    else:
        label.set_color('black')  # fallback color if missing


ax.set_xlabel("Metabolic Functions", fontsize=14)

pos = ax.get_position()


plt.tight_layout(rect=[0, 0, 0.9, 1])  # leave space on the right for colorbars
plt.show()
