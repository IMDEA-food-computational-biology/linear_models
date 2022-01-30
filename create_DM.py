from collections import Counter
import pandas as pd
from pyitlib import discrete_random_variable as drv
import numpy as np

def clean_DM(DM):

    """
    Cleans DM for use in linear models
    Removes columns with only one value
    Removes columns where each treatment has only one value
    """
    non_unique_cols = [c for c
        in list(DM)
        if len(DM[c].unique()) > 1]

    DM = DM[non_unique_cols]
    for var in DM.columns:

        if var == "treatment" or var == "GSM" or var == "color":

            continue

        if drv.entropy_conditional(DM[var].astype(str), DM["treatment"].astype(str)) == 0:

            DM = DM.drop(var, axis = 1)

    return DM

def one_hot_2_color(DM):

    """
    one hot encode the DM for two color arrays
    since our R scripts cant deal with those well
    """
    raise NotImplementedError("2 color studies still not implemented :(")

def write_to_file(accession, cursor, path = "./"):

    study_type_query = f"select study_type from study where accession = {accession}"
    cursor.execute(study_type_query)
    study_type = cursor.fetchall()[0][0]
    two_color = False
    if study_type == 'Two channel array':

        two_color = True
        query = f"""select sample.GSM, sample.treatment, sample.origin_name,
        sample.time_point, sample.concentration, misc_sample.attribute_value as 'color'
        from study natural join sample natural join misc_sample
        where study.accession = {accession}
        and attribute_name = 'color' """

    else:

        two_color = False
        query = f"""select sample.GSM, sample.treatment, sample.origin_name, sample.time_point, sample.concentration
        from sample natural join study where study.accession =  {accession}"""

    cursor.execute(query)
    result = cursor.fetchall()
    column_names = cursor.column_names
    DM = pd.DataFrame(result, columns=column_names)
    DM = clean_DM(DM)
    if two_color:
        
        check_2_per_GSM = set(Counter(DM.GSM).values())
        if len(check_2_per_GSM) != 1:

            assert False

        if 2 not in check_2_per_GSM:

            assert False
        
        #remove reference samples, as they are always redundant
        
        #check if the experiment had dye swaps
        DM = DM[DM["treatment"] != "ref"]
        colors = set(DM["color"])

        if len(set(DM["color"])) == 1:

            if colors == {"Cy5"}:
                DM = DM.drop(["color"],1)
                #If no dye swap study essentially boils down to a one color study
                study_type = "One color array"
                
            else:
                
                raise NotImplementedError("No implementation for 2 color refernce studies with all being Cy3")
       

    with open( path + str(accession) + "_DM.tsv", "w+") as f:

        f.write("#" + study_type + "\n")
        DM.to_csv(f, sep = "\t", index = 0)

    node_query = f"select sample.GSM, nodes.node_id from study natural join sample natural join nodes where study.accession = {accession}"
    cursor.execute(node_query)
    result = cursor.fetchall()
    column_names = cursor.column_names
    node_ids = pd.DataFrame(result, columns=column_names)
    node_ids.to_csv(path + accession + "_" + "node_ids.tsv", sep = "\t", index = 0)

if __name__ == "__main__":

    from getpass import getpass
    from mysql.connector import connect, Error
    import sys

    gses = sys.argv[1]
    try:

        connection = connect(host = "localhost",
                     user = input("Enter username: "),
                    password = getpass("Enter password: "),
                    database = "FooDrugs")
        cursor = connection.cursor()
        with open(gses) as f:

            for line in f:

                print(line[:-1])
                gse = line[:-1]
                write_to_file(gse,cursor, path = "./DM/")

    finally:

        connection.close()
