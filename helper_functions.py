import csv, xlwt

def integrate_up_down_csv_files_into_excel(outfilename):
    wb = xlwt.Workbook()

    filenames = ['up_up.csv', 'up_down.csv', 'down_up.csv', 'down_down.csv']

    for filename in filenames:
        f_short_name = filename.split(".csv")[0]
        ws = wb.add_sheet(f_short_name.title())
        data_for_excel = csv.reader(open(filename, "r"))
        for rowx, row in enumerate(data_for_excel):
            for colx, value in enumerate(row):
                ws.write(rowx, colx, value)
                
    wb.save(outfilename)


def parse_diff_list_table(filename, which_p_value):
    """
    This script is for parsing 
    diff_list_table.tsv
    to get a better understanding 
    of interaction term and 
    multifactor design analysis.
    """
    type_of_genes = {"up_up": [], "down_down": [], "down_up": [], "up_down": []}

    header = " "

    # Get header
    with open(filename, "r") as f:
        for lines in f:
            lines = lines.strip().split("\t")
            header = lines
            break

    with open(filename, "r") as f:
        next(f)
        for lines in f:
            lines = lines.strip().split("\t")
            gene = lines[0].split(",")[0]

            logfc1 = 0 if lines[1] == "NA" else float(lines[1])
            logfc2 = 0 if lines[4]  == "NA" else float(lines[4])
            interaction_fc = 0 if lines[7] == "NA" else float(lines[7])

            if which_p_value == "pvalue":
                pvalue1 = 1 if lines[2] == "NA" else float(lines[2])
                pvalue2 = 1 if lines[5] == "NA" else float(lines[5])
                interaction_pvalue = 1 if lines[5] == "NA" else float(lines[8])
                newfile_header = ",".join(["GeneName", header[0], header[1], header[3], \
                                            header[4], header[6], header[7]])

            elif which_p_value == "padj":
                pvalue1 = 1 if lines[3] == "NA" else float(lines[3])
                pvalue2 = 1 if lines[6] == "NA" else float(lines[6])
                interaction_pvalue = 1 if lines[9] == "NA" else float(lines[9])
                newfile_header = ",".join(["GeneName", header[0], header[2], header[3], \
                                            header[5], header[6], header[8]])   

            if interaction_pvalue < 0.05:
                list_of_values = [gene, logfc1, pvalue1, logfc2, pvalue2, interaction_fc, interaction_pvalue]

                if logfc1 >= 0.0 and interaction_fc >= 0.0:
                    type_of_genes["up_up"].append(map(str, list_of_values))

                elif logfc1 < 0.0 and interaction_fc < 0.0:
                    type_of_genes["down_down"].append(map(str, list_of_values))

                elif logfc1 < 0.0 and interaction_fc >= 0.0:
                    type_of_genes["down_up"].append(map(str, list_of_values))

                elif logfc1 >= 0.0 and interaction_fc < 0.0:
                    type_of_genes["up_down"].append(map(str, list_of_values))

    for key in type_of_genes:
        with open(key+".csv", "w") as outfile:
            outfile.write(newfile_header+"\n")
            for gene in type_of_genes[key]: 
                outfile.write(",".join(gene)+"\n")

    return(type_of_genes)
