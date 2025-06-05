import csv
import numpy as np

# File path
#file_path = 'R:/Gre25/Summary/thicknesses/thicknesses_vs_time.csv'


def open_csv_table(file_path):
    # Flag to skip the first 2 rows (the first one are the caterories and the second one are the types of each column)
    skipfirst2rows = True

    # Dictionaries to store categories and lists
    list_categories = []
    list_types = []
    dict_lists = {}

    # Open and read the CSV file
    with open(file_path, 'r') as file:
        reader = csv.reader(file,delimiter=';')
        count = 0

        for row in reader:
            if count == 0 and skipfirst2rows:  # Handle the header row
                list_categories = [row[i] for i in range(len(row))]  # Map index to column name
                for i in range(len(row)):
                    dict_lists[row[i]] = []  # Initialize empty lists for each column
                count += 1
            elif count == 1 and skipfirst2rows:
                list_types = [row[i] for i in range(len(row))]  # Map index to column name
                count+=1
            else:
                for i in range(len(row)):
                    if list_types[i]=='int':
                        dict_lists[list_categories[i]].append(int(row[i]))  # Append each cell to its corresponding list
                    elif list_types[i]=='float':
                        dict_lists[list_categories[i]].append(float(row[i].replace(',', '.')))
                    else: # par defaut ca sera str
                        dict_lists[list_categories[i]].append(row[i])
                count += 1
    
    return (dict_lists , list_categories , list_types)

# Output the results
#print("Categories (column headers):", list_categories)
#print("Categories (column headers):", list_types)
#print("Data (lists for each column):",dict_lists)






# deuxieme variante de la fonction, plus adaptée au cas de la manip de grenoble

def open_csv_table_experiments(file_path):
    # Flag to skip the first 2 rows (the first one are the caterories and the second one are the types of each column)
    skipfirst2rows = True

    # Dictionaries to store categories and lists
    list_categories = []
    list_types = []
    dict_lists = {}

    # Open and read the CSV file
    with open(file_path, 'r') as file:
        reader = csv.reader(file,delimiter=';')
        count = 0

        for row in reader:
            if count == 0 and skipfirst2rows:  # Handle the header row
                list_categories = [row[i] for i in range(len(row))]  # Map index to column name
                for i in range(len(row)):
                    dict_lists[row[i]] = []  # Initialize empty lists for each column
                count += 1
            elif count == 1 and skipfirst2rows:
                list_types = [row[i] for i in range(len(row))]  # Map index to column name
                count+=1
            else:
                for i in range(len(row)):
                    if list_types[i]=='int':
                        if (row[i]=='nan') or (row[i]==''):
                            dict_lists[list_categories[i]].append(np.nan)
                        else:
                            dict_lists[list_categories[i]].append(int(row[i]))  # Append each cell to its corresponding list
                    elif list_types[i]=='float':
                        if (row[i]=='nan') or (row[i]==''):
                            dict_lists[list_categories[i]].append(np.nan)
                        else:
                            dict_lists[list_categories[i]].append(float(row[i].replace(',', '.')))
                    else: # par defaut ca sera str
                        dict_lists[list_categories[i]].append(row[i])
                count += 1
    
    return (dict_lists , list_categories , list_types)
