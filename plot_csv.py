import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import input_data

csv_directory = "."  # Replace with the appropriate directory if not the current one

# Find all CSV files ending with '_properties.csv'
csv_files = [file for file in os.listdir(csv_directory) if file.endswith('_properties.csv')]

if not csv_files:
    print("No CSV files ending with '_properties.csv' found.")
else:
    for csv_file in csv_files:
        plt.plot(); plot = True
        #extract location:
        location = csv_file[ csv_file.find('-'):  csv_file.find('-',  csv_file.find('-') +1 ) -1 ]
        print('Location: ', location)

        # Read the CSV file into a DataFrame
        file_path = os.path.join(csv_directory, csv_file)
        df = pd.read_csv(file_path)

        # Ensure the required columns are present
        required_columns = {'Y', 'Z'}
        if not required_columns.issubset(df.columns):
            print(f"Skipping {csv_file}: Missing required columns {required_columns}.")
            continue

        # Calculate the new column R
        df['R'] = np.sqrt(df['Y']**2 + df['Z']**2)

        # Sort the DataFrame by R
        df = df.sort_values(by='R').reset_index(drop=True)

        # Drop excluded columns from the columns to plot
        excluded_columns = {'X', 'Y', 'Z', 'breakup', 'coalesce', 'R'}
        columns_to_plot = [col for col in df.columns if col not in excluded_columns]

        if not columns_to_plot:
            print(f"Skipping {csv_file}: No columns available for plotting against R.")
            continue
        
        # Plot each column versus R
        plt.figure()
        for column in columns_to_plot:         

            if 'r_0' in column:
                plt.figure()
                plt.scatter(df['R']*1000, df[column], label=r'$r_0$')
                plot = False
            elif 'L' == column:
                plt.scatter(df['R']*1000, df[column], label=r'$L$', marker='^')
                plot = False
            elif 'B' == column:
                plt.scatter(df['R']*1000, df[column], label=r'$B$', marker='*')
                plt.ylabel("Ellipsoid")
                plot = True
            elif "fibrillation" in column:
                plt.figure()
                plt.scatter(df['R']*1000, df[column], label=r'$FI$')
                plt.ylabel(r'$FI$: Fibrillation index')
                plt.title(r'$FI$ vs $R$')
                plot = True

            elif "N_11" in column:
                plt.figure()
                plt.scatter(df['R']*1000, df[column], label=r'$N_{11}$', marker='s')
                plt.ylabel(r'Stress $[Pa]$')
                plt.title(r'Stress vs $R$')
                plot = False
            elif "shear" in column:
                plt.scatter(df['R']*1000, df[column], label=r'$\tau_{xr}$')
                plot = True

            elif "Ca_e" == column:
                plt.figure()
                plt.scatter(df['R']*1000, df[column], label=r'$Ca_e$', marker='s')
                plt.ylabel(r'$Ca$')
                plt.title(r'$Ca$ vs $R$')
                plot = False
            elif "Ca_s" == column:
                plt.scatter(df['R']*1000, df[column], label=r'$Ca_s$')
                plot = True

            elif "Ca_star_e" in column:
                plt.figure()
                plt.scatter(df['R']*1000, df[column], label=r'$Ca_e*$', marker='s')
                plt.title(r'$Ca^*$ vs $R$')
                plt.ylabel(r'$Ca*$')
                plot = False
            elif "Ca_star_s" in column:
                plt.scatter(df['R']*1000, df[column], label=r'$Ca_s*$')
                plot = True

            else: 
                plt.figure()
                plt.scatter(df['R']*1000, df[column], label=column)
                plt.ylabel(column)
                plt.title(column + r' vs $R$')
                plot = True

            plt.xlabel(r'$R[mm]$')
            plt.legend()
            plt.grid()
            if plot: plt.savefig(input_data.filename + column + '_' + location+'.png', dpi=300)
            plot = True


print("All files processed.")